/* RODAS3 for a small dense system (protolife/rodas.h), after Genesis core/rosenbrock.c and core/dense.c. */
#include "protolife/rodas.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "diagnostic_internal.h"

/* KPP's form: G = I/(h gamma) - J; G K_i = f(Y_i) + sum_j C_ij K_j / h; Y_i = y + sum_j A_ij K_j;
 * y' = y + 2 K1 + K3 + K4; the error estimate is K4. */
#define GAMMA 0.5
static const double C[4][3] = {{0, 0, 0}, {4.0, 0, 0}, {1.0, -1.0, 0}, {1.0, -1.0, -8.0 / 3.0}};
static const double A[4][3] = {{0, 0, 0}, {0, 0, 0}, {2.0, 0, 0}, {2.0, 0, 1.0}};

/* LU with partial pivoting of a row-major n x n matrix, in place; order records the row swaps. -1 if singular. */
static int lu_factorise(int n, double *a, int *order) {
    for (int c = 0; c < n; c++) {
        int best = c;
        for (int r = c + 1; r < n; r++) best = fabs(a[r * n + c]) > fabs(a[best * n + c]) ? r : best;
        if (a[best * n + c] == 0.0) return -1;
        order[c] = best;
        for (int k = 0; k < n; k++) {
            double swap = a[c * n + k];
            a[c * n + k] = a[best * n + k], a[best * n + k] = swap;
        }
        for (int r = c + 1; r < n; r++) {
            a[r * n + c] /= a[c * n + c];
            for (int k = c + 1; k < n; k++) a[r * n + k] -= a[r * n + c] * a[c * n + k];
        }
    }
    return 0;
}

/* Solves a x = b with the factorised a; b is replaced by x. Every row swap first: the factorisation swapped whole
 * rows, multipliers included, so L is stored in the final order of the rows. */
static void lu_solve(int n, const double *a, const int *order, double *b) {
    for (int c = 0; c < n; c++) {
        double swap = b[c];
        b[c] = b[order[c]], b[order[c]] = swap;
    }
    for (int c = 0; c < n; c++)
        for (int r = c + 1; r < n; r++) b[r] -= a[r * n + c] * b[c];
    for (int r = n - 1; r >= 0; r--) {
        for (int k = r + 1; k < n; k++) b[r] -= a[r * n + k] * b[k];
        b[r] /= a[r * n + r];
    }
}

typedef struct {
    int *order;
    double *jac, *g, *k[4], *stage, *f, *state, *next;
    int callback;  /* the last callback's nonzero return, 0 when none */
} work;

static int work_alloc(work *w, int n) {
    size_t square = (size_t)n * (size_t)n;
    w->order = malloc(sizeof *w->order * (size_t)n);
    w->jac = malloc(sizeof(double) * (2 * square + 8 * (size_t)n));
    if (w->order == NULL || w->jac == NULL) return free(w->order), free(w->jac), -1;
    w->g = w->jac + square;
    double *rest = w->g + square;
    for (int s = 0; s < 4; s++) w->k[s] = rest + (size_t)s * (size_t)n;
    w->stage = rest + 4 * (size_t)n, w->f = rest + 5 * (size_t)n, w->state = rest + 6 * (size_t)n;
    w->next = rest + 7 * (size_t)n;
    w->callback = 0;
    return 0;
}

/* One trial step of length h from y into next: the scaled error (accept if <= 1), INFINITY for a refused state. */
static double trial(const pl_rodas_options *o, work *w, double t, const double *y, double h, double *next) {
    int n = o->size, *order = w->order;
    double *jac = w->jac, *g = w->g, *stage = w->stage, *f = w->f, **k = w->k;
    if ((w->callback = o->jacobian(t, y, jac, o->context)) != 0) return INFINITY;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) g[i * n + j] = (i == j) / (h * GAMMA) - jac[i * n + j];
    if (lu_factorise(n, g, order) != 0) return INFINITY;
    for (int s = 0; s < 4; s++) {
        if (s != 1) {  /* stage 2 reuses stage 1's derivative (A21 = 0) */
            for (int i = 0; i < n; i++) {
                stage[i] = y[i];
                for (int j = 0; j < s; j++) stage[i] += A[s][j] * k[j][i];
            }
            if ((w->callback = o->rhs(t, stage, f, o->context)) != 0) return INFINITY;
        }
        for (int i = 0; i < n; i++) {
            k[s][i] = f[i];
            for (int j = 0; j < s; j++) k[s][i] += C[s][j] * k[j][i] / h;
        }
        lu_solve(n, g, order, k[s]);
    }
    double worst = 0.0;
    for (int i = 0; i < n; i++) {
        next[i] = y[i] + 2.0 * k[0][i] + k[2][i] + k[3][i];
        double allowance = o->tolerance * fmax(fmax(fabs(next[i]), fabs(y[i])), o->floor[i]);
        if (!isfinite(next[i]) || next[i] < -o->floor[i]) return INFINITY;  /* below zero by more than the noise */
        worst = fmax(worst, fabs(k[3][i]) / allowance);
    }
    return worst;
}

static pl_status check(const pl_rodas_options *o, const double *y, double duration, const double *step,
                       pl_diagnostic *diagnostic) {
    if (o == NULL || y == NULL || step == NULL || o->size < 1 || o->rhs == NULL || o->jacobian == NULL ||
        o->floor == NULL)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "RODAS3 needs a size, f, its Jacobian, floors, y and a step");
    if (!(o->tolerance > 0.0) || !(o->smallest > 0.0) || o->max_trials < 0 || !(duration >= 0.0) || !isfinite(duration))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "RODAS3 needs a tolerance, a smallest step and a duration > 0");
    for (int i = 0; i < o->size; i++)
        if (!(o->floor[i] > 0.0) || !isfinite(y[i]))
            return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "component %d needs a floor > 0 and a finite value", i);
    return PL_OK;
}

pl_status pl_rodas_advance(const pl_rodas_options *o, double t, double *y, double duration, double *step_inout,
                           pl_ode_statistics *statistics, pl_diagnostic *diagnostic) {
    pl_status status = check(o, y, duration, step_inout, diagnostic);
    if (status != PL_OK) return status;
    work w;
    if (work_alloc(&w, o->size) != 0) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, t, "no memory for RODAS3");
    memcpy(w.state, y, sizeof(double) * (size_t)o->size);
    double time = 0.0, step = *step_inout > 0.0 ? *step_inout : duration, error = 0.0, h = step;
    long trials = 0, accepted = 0, refused = 0;
    const char *why = NULL;
    while (time < duration && why == NULL) {
        if (o->max_trials > 0 && trials++ >= o->max_trials) {
            why = "it took the most trial steps allowed";
            break;
        }
        h = fmin(step, duration - time);
        error = trial(o, &w, t + time, w.state, h, w.next);
        if (error <= 1.0) {
            memcpy(w.state, w.next, sizeof(double) * (size_t)o->size);
            time = h < duration - time ? time + h : duration;
            accepted++;
        } else {
            refused++;
        }
        double factor = !isfinite(error) ? 0.25 : error > 0.0 ? 0.9 * pow(error, -1.0 / 3.0) : 5.0;
        step = h * fmin(fmax(factor, 0.2), 5.0);
        if (!(error <= 1.0) && !(step >= o->smallest)) why = "its step fell below the smallest allowed";
    }
    if (statistics != NULL) statistics->steps += accepted, statistics->error_test_failures += refused;
    if (why == NULL) {
        memcpy(y, w.state, sizeof(double) * (size_t)o->size);
        *step_inout = step;
        pl_succeed(diagnostic);
    } else {
        status = w.callback != 0 ? PL_ERROR_MODEL : PL_ERROR_SOLVER;
        pl_fail(diagnostic, status, w.callback, t + time,
                "RODAS3 stopped at t = %.17g on the way from %.17g to %.17g: %s (last trial step %.3g s, scaled error "
                "%.3g%s)", t + time, t, t + duration, why, h, error, w.callback != 0 ? "; a callback refused" : "");
    }
    free(w.order), free(w.jac);
    return status;
}
