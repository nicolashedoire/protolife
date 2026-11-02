/* Robertson's stiff chemistry (1966) with pl_ode, printed as CSV on a logarithmic time grid: time (s), A, B, C, and the
 * integrator's cumulative step count. Three reactions whose rate constants span nine orders of magnitude.
 *
 *   robertson [last_time_s] [points_per_decade]
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "protolife/protolife.h"

static int rhs(double t, const double *y, double *dydt, void *context) {
    (void)t, (void)context;
    dydt[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
    dydt[2] = 3e7 * y[1] * y[1];
    dydt[1] = -dydt[0] - dydt[2];
    return 0;
}

static int jacobian(double t, const double *y, double *j, void *context) {
    (void)t, (void)context;
    j[0] = -0.04, j[1] = 1e4 * y[2], j[2] = 1e4 * y[1];
    j[3] = 0.04, j[4] = -1e4 * y[2] - 6e7 * y[1], j[5] = -1e4 * y[1];
    j[6] = 0.0, j[7] = 6e7 * y[1], j[8] = 0.0;
    return 0;
}

int main(int argc, char **argv) {
    double last = argc > 1 ? atof(argv[1]) : 4e10;
    int per_decade = argc > 2 ? atoi(argv[2]) : 20;
    double y0[3] = {1.0, 0.0, 0.0}, absolute[3] = {1e-14, 1e-18, 1e-14};
    pl_ode_options o = {.size = 3, .rhs = rhs, .jacobian = jacobian, .relative_tolerance = 1e-11,
                        .absolute_tolerance = absolute, .nonnegative = 1, .max_steps = 100000};
    pl_ode *ode = NULL;
    pl_diagnostic d;
    if (pl_ode_create(&o, 0.0, y0, &ode, &d) != PL_OK) return fprintf(stderr, "%s\n", d.message), 1;
    printf("time,A,B,C,steps\n");
    for (int k = 0;; k++) {
        double t = 1e-6 * pow(10.0, (double)k / per_decade);
        if (t > last * (1.0 + 1e-12)) break;
        if (pl_ode_advance(ode, t, NULL, &d) != PL_OK) {
            fprintf(stderr, "%s: %s\n", pl_status_name(d.status), d.message);
            pl_ode_free(ode);
            return 1;
        }
        pl_ode_statistics s;
        pl_ode_get_statistics(ode, &s);
        const double *y = pl_ode_state(ode);
        printf("%.6g,%.9g,%.9g,%.9g,%ld\n", t, y[0], y[1], y[2], s.steps);
    }
    pl_ode_free(ode);
    return 0;
}
