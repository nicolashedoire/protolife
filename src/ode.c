/* Stiff ODE integration over CVODE (protolife/ode.h). */
#include "protolife/ode.h"

#include <cvode/cvode.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "diagnostic_internal.h"

struct pl_ode {
    pl_ode_options options;    /* absolute_tolerance points to `absolute`, owned here */
    double *absolute;
    double *jacobian;          /* the row-major work array handed to the Jacobian callback; NULL without one */
    double t;                  /* the time of `y` */
    SUNContext context;
    void *cvode;
    N_Vector y;                /* the state the caller sees */
    N_Vector work;             /* CVODE's output during a call, copied to y on success */
    N_Vector tolerances, constraints;
    SUNMatrix matrix;
    SUNLinearSolver solver;
    pl_ode_statistics before;  /* counts of the runs that ended with a reinitialisation */
    char note[256];            /* the last callback failure */
    char sundials[256];        /* the last message of SUNDIALS's error handler */
};

static int first_not_finite(const double *values, int n) {
    for (int i = 0; i < n; i++)
        if (!isfinite(values[i])) return i;
    return -1;
}

static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void *data) {
    pl_ode *ode = data;
    double *dydt = N_VGetArrayPointer(ydot);
    int result = ode->options.rhs(t, N_VGetArrayPointer(y), dydt, ode->options.context);
    if (result != 0) {
        snprintf(ode->note, sizeof ode->note, "the right-hand side returned %d at t = %.17g", result, t);
        return result;
    }
    int bad = first_not_finite(dydt, ode->options.size);
    if (bad >= 0) {
        snprintf(ode->note, sizeof ode->note, "the right-hand side is not finite for component %d at t = %.17g", bad,
                 t);
        return 1;
    }
    return 0;
}

static int jacobian(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix matrix, void *data, N_Vector work1,
                    N_Vector work2, N_Vector work3) {
    (void)fy, (void)work1, (void)work2, (void)work3;
    pl_ode *ode = data;
    int n = ode->options.size;
    int result = ode->options.jacobian(t, N_VGetArrayPointer(y), ode->jacobian, ode->options.context);
    if (result != 0) {
        snprintf(ode->note, sizeof ode->note, "the Jacobian returned %d at t = %.17g", result, t);
        return result;
    }
    int bad = first_not_finite(ode->jacobian, n * n);
    if (bad >= 0) {
        snprintf(ode->note, sizeof ode->note, "the Jacobian is not finite at row %d, column %d, t = %.17g", bad / n,
                 bad % n, t);
        return 1;
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) SM_ELEMENT_D(matrix, i, j) = ode->jacobian[i * n + j];
    return 0;
}

static int roots(sunrealtype t, N_Vector y, sunrealtype *g, void *data) {
    pl_ode *ode = data;
    int result = ode->options.event(t, N_VGetArrayPointer(y), g, ode->options.context);
    if (result != 0) {
        snprintf(ode->note, sizeof ode->note, "the event functions returned %d at t = %.17g", result, t);
        return result;
    }
    int bad = first_not_finite(g, ode->options.events);
    if (bad >= 0) {
        snprintf(ode->note, sizeof ode->note, "event function %d is not finite at t = %.17g", bad, t);
        return -1;
    }
    return 0;
}

static void keep_message(int line, const char *function, const char *file, const char *message, SUNErrCode code,
                         void *data, SUNContext context) {
    (void)line, (void)function, (void)file, (void)code, (void)context;
    pl_ode *ode = data;
    snprintf(ode->sundials, sizeof ode->sundials, "%s", message != NULL ? message : "");
}

static void add_counts(const pl_ode *ode, pl_ode_statistics *s) {
    long steps = 0, rhs_evaluations = 0, linear_rhs = 0, jacobians = 0, error_tests = 0, convergence = 0;
    CVodeGetNumSteps(ode->cvode, &steps);
    CVodeGetNumRhsEvals(ode->cvode, &rhs_evaluations);
    CVodeGetNumLinRhsEvals(ode->cvode, &linear_rhs);
    CVodeGetNumJacEvals(ode->cvode, &jacobians);
    CVodeGetNumErrTestFails(ode->cvode, &error_tests);
    CVodeGetNumNonlinSolvConvFails(ode->cvode, &convergence);
    s->steps += steps, s->rhs_evaluations += rhs_evaluations + linear_rhs, s->jacobian_evaluations += jacobians;
    s->error_test_failures += error_tests, s->convergence_failures += convergence;
}

/* Reinitialises CVODE at (t, y), keeping the counts of the run it ends; the first step tried is the last one taken. */
static int reinitialise(pl_ode *ode, double t, N_Vector y) {
    double last = 0.0;
    CVodeGetLastStep(ode->cvode, &last);
    add_counts(ode, &ode->before);
    int flag = CVodeReInit(ode->cvode, t, y);
    if (flag == CV_SUCCESS && last > 0.0) flag = CVodeSetInitStep(ode->cvode, last);
    return flag;
}

static pl_status cvode_failure(pl_ode *ode, int flag, int linear, const char *what, double t,
                               pl_diagnostic *diagnostic) {
    char *name = linear ? CVodeGetLinReturnFlagName(flag) : CVodeGetReturnFlagName(flag);
    pl_status status = pl_fail(diagnostic, flag == CV_MEM_FAIL ? PL_ERROR_MEMORY : PL_ERROR_ARGUMENT, flag, t,
                               "%s failed with %s: %s", what, name != NULL ? name : "an unknown flag", ode->sundials);
    free(name);
    return status;
}

void pl_ode_free(pl_ode *ode) {
    if (ode == NULL) return;
    if (ode->cvode != NULL) CVodeFree(&ode->cvode);
    if (ode->solver != NULL) SUNLinSolFree(ode->solver);
    if (ode->matrix != NULL) SUNMatDestroy(ode->matrix);
    N_Vector vectors[] = {ode->y, ode->work, ode->tolerances, ode->constraints};
    for (size_t i = 0; i < sizeof vectors / sizeof *vectors; i++)
        if (vectors[i] != NULL) N_VDestroy(vectors[i]);
    if (ode->context != NULL) SUNContext_Free(&ode->context);
    free(ode->absolute);
    free(ode->jacobian);
    free(ode);
}

static pl_status check_options(const pl_ode_options *o, double t0, const double *y0, pl_diagnostic *diagnostic) {
    if (o == NULL || y0 == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the options and y0 are required");
    if (o->size <= 0) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the size must be positive, not %d", o->size);
    if (o->rhs == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the right-hand side is required");
    if (!(o->relative_tolerance > 0.0) || !isfinite(o->relative_tolerance))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the relative tolerance must be positive and finite");
    if (o->absolute_tolerance == NULL)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the absolute tolerances are required");
    for (int i = 0; i < o->size; i++)
        if (!(o->absolute_tolerance[i] > 0.0) || !isfinite(o->absolute_tolerance[i]))
            return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN,
                           "the absolute tolerance of component %d must be positive and finite", i);
    if (o->max_steps < 0 || !(o->min_step >= 0.0) || !(o->max_step >= 0.0) || !(o->initial_step >= 0.0) ||
        o->max_nonlinear_iterations < 0 || o->max_convergence_failures < 0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the step limits must not be negative");
    if (o->events < 0 || (o->events > 0 && o->event == NULL))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "events > 0 needs the event functions");
    if (!isfinite(t0)) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the initial time must be finite");
    int bad = first_not_finite(y0, o->size);
    if (bad >= 0) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, t0, "component %d of y0 is not finite", bad);
    if (o->nonnegative)
        for (int i = 0; i < o->size; i++)
            if (y0[i] < 0.0)
                return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, t0, "component %d of y0 is negative (%.17g)", i, y0[i]);
    return PL_OK;
}

/* Every CVODE setting of a new integrator; on failure the caller frees what was made. */
static pl_status configure(pl_ode *ode, double t0, pl_diagnostic *diagnostic) {
    const pl_ode_options *o = &ode->options;
    ode->cvode = CVodeCreate(CV_BDF, ode->context);
    ode->matrix = SUNDenseMatrix(o->size, o->size, ode->context);
    ode->solver = ode->matrix != NULL ? SUNLinSol_Dense(ode->y, ode->matrix, ode->context) : NULL;
    if (ode->cvode == NULL || ode->matrix == NULL || ode->solver == NULL)
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, t0, "CVODE could not be allocated: %s", ode->sundials);
    const char *what = "CVodeInit";
    int linear = 0, flag = CVodeInit(ode->cvode, rhs, t0, ode->y);
    if (flag == CV_SUCCESS) what = "CVodeSetUserData", flag = CVodeSetUserData(ode->cvode, ode);
    if (flag == CV_SUCCESS)
        what = "CVodeSVtolerances", flag = CVodeSVtolerances(ode->cvode, o->relative_tolerance, ode->tolerances);
    if (flag == CV_SUCCESS)
        what = "CVodeSetLinearSolver", linear = 1, flag = CVodeSetLinearSolver(ode->cvode, ode->solver, ode->matrix);
    if (flag == CV_SUCCESS && o->jacobian != NULL) what = "CVodeSetJacFn", flag = CVodeSetJacFn(ode->cvode, jacobian);
    if (flag == CV_SUCCESS) linear = 0;
    if (flag == CV_SUCCESS && o->nonnegative)
        what = "CVodeSetConstraints", flag = CVodeSetConstraints(ode->cvode, ode->constraints);
    if (flag == CV_SUCCESS && o->max_steps > 0)
        what = "CVodeSetMaxNumSteps", flag = CVodeSetMaxNumSteps(ode->cvode, o->max_steps);
    if (flag == CV_SUCCESS && o->min_step > 0.0) what = "CVodeSetMinStep", flag = CVodeSetMinStep(ode->cvode, o->min_step);
    if (flag == CV_SUCCESS && o->max_step > 0.0) what = "CVodeSetMaxStep", flag = CVodeSetMaxStep(ode->cvode, o->max_step);
    if (flag == CV_SUCCESS && o->max_nonlinear_iterations > 0)
        what = "CVodeSetMaxNonlinIters", flag = CVodeSetMaxNonlinIters(ode->cvode, o->max_nonlinear_iterations);
    if (flag == CV_SUCCESS && o->max_convergence_failures > 0)
        what = "CVodeSetMaxConvFails", flag = CVodeSetMaxConvFails(ode->cvode, o->max_convergence_failures);
    if (flag == CV_SUCCESS && o->initial_step > 0.0)
        what = "CVodeSetInitStep", flag = CVodeSetInitStep(ode->cvode, o->initial_step);
    if (flag == CV_SUCCESS && o->events > 0) what = "CVodeRootInit", flag = CVodeRootInit(ode->cvode, o->events, roots);
    return flag == CV_SUCCESS ? PL_OK : cvode_failure(ode, flag, linear, what, t0, diagnostic);
}

pl_status pl_ode_create(const pl_ode_options *options, double t0, const double *y0, pl_ode **out,
                        pl_diagnostic *diagnostic) {
    if (out == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the integrator's address is required");
    *out = NULL;
    pl_status status = check_options(options, t0, y0, diagnostic);
    if (status != PL_OK) return status;
    int n = options->size;
    pl_ode *ode = calloc(1, sizeof *ode);
    if (ode == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, t0, "the integrator could not be allocated");
    ode->options = *options, ode->t = t0;
    ode->absolute = malloc((size_t)n * sizeof *ode->absolute);
    ode->jacobian = options->jacobian != NULL ? malloc((size_t)n * (size_t)n * sizeof *ode->jacobian) : NULL;
    if (ode->absolute == NULL || (options->jacobian != NULL && ode->jacobian == NULL) ||
        SUNContext_Create(SUN_COMM_NULL, &ode->context) != SUN_SUCCESS) {
        pl_ode_free(ode);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, t0, "the integrator could not be allocated");
    }
    memcpy(ode->absolute, options->absolute_tolerance, (size_t)n * sizeof *ode->absolute);
    ode->options.absolute_tolerance = ode->absolute;
    SUNContext_ClearErrHandlers(ode->context);
    SUNContext_PushErrHandler(ode->context, keep_message, ode);
    ode->y = N_VNew_Serial(n, ode->context);
    ode->work = ode->y != NULL ? N_VClone(ode->y) : NULL;
    ode->tolerances = ode->y != NULL ? N_VClone(ode->y) : NULL;
    ode->constraints = ode->y != NULL && options->nonnegative ? N_VClone(ode->y) : NULL;
    if (ode->work == NULL || ode->tolerances == NULL || (options->nonnegative && ode->constraints == NULL)) {
        pl_ode_free(ode);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, t0, "the integrator's vectors could not be allocated");
    }
    memcpy(N_VGetArrayPointer(ode->y), y0, (size_t)n * sizeof *y0);
    memcpy(N_VGetArrayPointer(ode->tolerances), ode->absolute, (size_t)n * sizeof *y0);
    if (ode->constraints != NULL) N_VConst(1.0, ode->constraints);  /* CVODE's code for "at least zero" */
    status = configure(ode, t0, diagnostic);
    if (status != PL_OK) {
        pl_ode_free(ode);
        return status;
    }
    *out = ode;
    pl_succeed(diagnostic);
    return PL_OK;
}

/* A failed call in which a callback refused is the model's failure, whatever flag CVODE ended with: refusals near the
 * edge of a model's domain usually make CVODE shrink its step until it runs out of steps. */
static pl_status failure_status(const pl_ode *ode, int flag) {
    if (flag == CV_MEM_FAIL) return PL_ERROR_MEMORY;
    if (ode->note[0] != '\0') return PL_ERROR_MODEL;
    switch (flag) {
    case CV_RHSFUNC_FAIL:
    case CV_FIRST_RHSFUNC_ERR:
    case CV_REPTD_RHSFUNC_ERR:
    case CV_UNREC_RHSFUNC_ERR:
    case CV_RTFUNC_FAIL: return PL_ERROR_MODEL;
    case CV_MEM_NULL:
    case CV_ILL_INPUT:
    case CV_BAD_T:
    case CV_TOO_CLOSE: return PL_ERROR_ARGUMENT;
    default: return PL_ERROR_SOLVER;
    }
}

pl_status pl_ode_advance(pl_ode *ode, double t_end, int *fired, pl_diagnostic *diagnostic) {
    if (ode == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the integrator is required");
    if (!isfinite(t_end) || !(t_end > ode->t))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, ode->t, "t_end (%.17g) must come after the current time",
                       t_end);
    if (ode->options.events > 0 && fired == NULL)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, ode->t, "an integrator with events needs `fired`");
    ode->note[0] = ode->sundials[0] = '\0';
    double reached = ode->t;
    int flag = CVodeSetStopTime(ode->cvode, t_end);
    if (flag == CV_SUCCESS) flag = CVode(ode->cvode, t_end, ode->work, &reached, CV_NORMAL);
    if (flag >= 0) {  /* CV_SUCCESS, CV_TSTOP_RETURN or CV_ROOT_RETURN */
        N_VScale(1.0, ode->work, ode->y);
        ode->t = reached;
        if (ode->options.events > 0) {
            if (flag == CV_ROOT_RETURN) CVodeGetRootInfo(ode->cvode, fired);
            else memset(fired, 0, (size_t)ode->options.events * sizeof *fired);
        }
        pl_succeed(diagnostic);
        return PL_OK;
    }
    double stopped = ode->t;
    CVodeGetCurrentTime(ode->cvode, &stopped);
    char *name = CVodeGetReturnFlagName(flag);
    pl_status status = failure_status(ode, flag);
    int restored = reinitialise(ode, ode->t, ode->y) == CV_SUCCESS;
    pl_fail(diagnostic, status, flag, stopped, "CVODE stopped with %s at t = %.17g on the way from %.17g to %.17g%s%s%s%s%s",
            name != NULL ? name : "an unknown flag", stopped, ode->t, t_end, ode->sundials[0] ? ": " : "", ode->sundials,
            ode->note[0] ? "; " : "", ode->note, restored ? "" : "; the integrator could not return to its last state");
    free(name);
    return status;
}

pl_status pl_ode_restart(pl_ode *ode, double t, const double *y, pl_diagnostic *diagnostic) {
    if (ode == NULL || y == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the integrator and y are required");
    pl_ode_options options = ode->options;
    pl_status status = check_options(&options, t, y, diagnostic);
    if (status != PL_OK) return status;
    ode->sundials[0] = '\0';
    memcpy(N_VGetArrayPointer(ode->work), y, (size_t)ode->options.size * sizeof *y);
    int flag = reinitialise(ode, t, ode->work);
    if (flag != CV_SUCCESS) {
        reinitialise(ode, ode->t, ode->y);
        return cvode_failure(ode, flag, 0, "CVodeReInit", t, diagnostic);
    }
    N_VScale(1.0, ode->work, ode->y);
    ode->t = t;
    pl_succeed(diagnostic);
    return PL_OK;
}

double pl_ode_time(const pl_ode *ode) { return ode != NULL ? ode->t : NAN; }

double pl_ode_last_step(const pl_ode *ode) {
    double last = 0.0;
    if (ode != NULL) CVodeGetLastStep(ode->cvode, &last);
    return last;
}

const double *pl_ode_state(const pl_ode *ode) { return ode != NULL ? N_VGetArrayPointer(ode->y) : NULL; }

int pl_ode_size(const pl_ode *ode) { return ode != NULL ? ode->options.size : 0; }

void pl_ode_get_statistics(const pl_ode *ode, pl_ode_statistics *statistics) {
    if (statistics == NULL) return;
    *statistics = (pl_ode_statistics){0};
    if (ode == NULL) return;
    *statistics = ode->before;
    add_counts(ode, statistics);
}
