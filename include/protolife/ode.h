#ifndef PROTOLIFE_ODE_H
#define PROTOLIFE_ODE_H

#include "protolife/diagnostic.h"

/*
 * Stiff ordinary differential equations dy/dt = f(t, y), integrated forward in time by CVODE (SUNDIALS: backward
 * differentiation formulas of variable order and step, Newton iterations, a dense linear solver), with optional event
 * functions whose zero crossings stop the integration.
 *
 * The integrator owns the state, and pl_ode_advance moves it forward only when the call succeeds. When the requested
 * time cannot be reached (too many steps, repeated error-test or convergence failures, a step below the minimum, a
 * component that would turn negative under `nonnegative`), or when a callback fails or returns a value that is not
 * finite, the call returns an error and the integrator goes back to the time and state it had before the call, ready
 * to be advanced again, possibly with other settings, or restarted. The status is PL_ERROR_MODEL when a callback
 * refused or returned a value that is not finite during the call, PL_ERROR_SOLVER otherwise; the diagnostic carries
 * CVODE's flag, the time CVODE reached, CVODE's own message and the callback's last failure. Nothing is clipped.
 */

typedef struct pl_ode pl_ode;

/* Every callback returns 0 on success; a positive value when the state lies outside the model's domain and a shorter
 * step may help (the integrator retries); a negative value for a failure no step can mend. */

/* f(t, y) into dydt (n values). */
typedef int (*pl_ode_rhs)(double t, const double *y, double *dydt, void *context);
/* The Jacobian, row-major n x n: jacobian[i * n + j] = d f_i / d y_j. */
typedef int (*pl_ode_jacobian)(double t, const double *y, double *jacobian, void *context);
/* The event functions g_k(t, y), k < events, into g: the integration stops at the first zero crossing of any. */
typedef int (*pl_ode_events)(double t, const double *y, double *g, void *context);

typedef struct pl_ode_options {
    int size;                          /* n, the number of components, > 0 */
    pl_ode_rhs rhs;                    /* required */
    pl_ode_jacobian jacobian;          /* NULL: CVODE approximates the Jacobian by differences */
    void *context;                     /* handed to every callback */
    double relative_tolerance;         /* > 0 */
    const double *absolute_tolerance;  /* n values > 0: per component, the size under which it is noise */
    int nonnegative;                   /* nonzero: a step that takes a component below zero is refused and retried
                                          shorter; the initial state must then be nonnegative */
    long max_steps;                    /* internal steps allowed per call of pl_ode_advance; 0 keeps CVODE's 500 */
    double min_step, max_step;         /* bounds on the internal step; 0 for none */
    int events;                        /* the number of event functions, 0 for none */
    pl_ode_events event;               /* required when events > 0 */
} pl_ode_options;

/* Counts since the integrator was created, over every call, successful or not. The right-hand side evaluations
 * include those made to approximate the Jacobian by differences. */
typedef struct pl_ode_statistics {
    long steps, rhs_evaluations, jacobian_evaluations, error_test_failures, convergence_failures;
} pl_ode_statistics;

/* An integrator at time t0 and state y0 (n values, copied). The options are copied too; the context is not. */
pl_status pl_ode_create(const pl_ode_options *options, double t0, const double *y0, pl_ode **ode,
                        pl_diagnostic *diagnostic);

/* Integrates up to t_end, later than the current time. On success the state is the one at the time reached, which is
 * t_end unless an event fired first; fired (`events` values, NULL when there are none) then holds +1 or -1 for each
 * event function that crossed zero upward or downward, and 0 for the others. */
pl_status pl_ode_advance(pl_ode *ode, double t_end, int *fired, pl_diagnostic *diagnostic);

/* Starts again from time t and state y (copied): after the caller changes the state between two calls, as a feeding,
 * a dilution or a division does. On failure the integrator keeps its previous time and state. */
pl_status pl_ode_restart(pl_ode *ode, double t, const double *y, pl_diagnostic *diagnostic);

double pl_ode_time(const pl_ode *ode);
/* The n values of the state, owned by the integrator and valid until it is freed. */
const double *pl_ode_state(const pl_ode *ode);
int pl_ode_size(const pl_ode *ode);
void pl_ode_get_statistics(const pl_ode *ode, pl_ode_statistics *statistics);
void pl_ode_free(pl_ode *ode);

#endif
