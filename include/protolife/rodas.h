#ifndef PROTOLIFE_RODAS_H
#define PROTOLIFE_RODAS_H

#include "protolife/diagnostic.h"
#include "protolife/ode.h"

/*
 * RODAS3, a Rosenbrock method for stiff autonomous systems dy/dt = f(y) (Sandu et al., Atmos. Environ. 31, 3459, 1997;
 * the coefficients of KPP's Rodas3): L-stable and stiffly accurate, order 3 with an embedded order-2 estimate, with the
 * Jacobian the caller gives and a dense LU solve. After Genesis core/rosenbrock.c.
 *
 * Unlike pl_ode's BDF, a Rosenbrock step solves linear systems only and iterates nothing: it suits models whose
 * right-hand side has kinks (a switch between two regimes), where a Newton iteration stalls. pl_vesicles uses it.
 *
 * A trial step is refused, and retried shorter, when its error is too large, when its state is not finite, or when a
 * component falls below minus its floor (the size under which it is noise): nothing is clipped. The error of component
 * i is measured against tolerance * max(|y_i|, floor_i). The state changes only when the whole duration is covered.
 */

typedef struct pl_rodas_options {
    int size;                   /* n, > 0 */
    pl_ode_rhs rhs;             /* f(y) into dydt; t is passed as the time of the step's start */
    pl_ode_jacobian jacobian;   /* required: row-major n x n */
    void *context;
    double tolerance;           /* relative error allowed per step, > 0 */
    const double *floor;        /* n values > 0: the size under which a component is noise */
    double smallest;            /* the shortest step allowed (> 0) */
    long max_trials;            /* trial steps one call may take; 0 for no limit */
} pl_rodas_options;

/* Integrates y (n values, from time t) over `duration` s. *step is the first trial step (0 starts from the whole
 * duration) and receives the next one. Counts of accepted and refused steps are added to statistics->steps and
 * statistics->error_test_failures when statistics is not NULL. On failure y is unchanged. */
pl_status pl_rodas_advance(const pl_rodas_options *options, double t, double *y, double duration, double *step,
                           pl_ode_statistics *statistics, pl_diagnostic *diagnostic);

#endif
