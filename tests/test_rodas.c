/* RODAS3 (protolife/rodas.h): an exact solution, Genesis's reference, and refusals that change nothing. */
#include <math.h>
#include <string.h>

#include "check.h"
#include "protolife/protolife.h"

static int decay(double t, const double *y, double *dydt, void *context) {
    (void)t, (void)context;
    dydt[0] = -2.0 * y[0];
    return 0;
}

static int decay_jacobian(double t, const double *y, double *j, void *context) {
    (void)t, (void)y, (void)context;
    j[0] = -2.0;
    return 0;
}

static int robertson(double t, const double *y, double *dydt, void *context) {
    (void)t, (void)context;
    dydt[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
    dydt[2] = 3e7 * y[1] * y[1];
    dydt[1] = -dydt[0] - dydt[2];
    return 0;
}

static int robertson_jacobian(double t, const double *y, double *j, void *context) {
    (void)t, (void)context;
    j[0] = -0.04, j[1] = 1e4 * y[2], j[2] = 1e4 * y[1];
    j[3] = 0.04, j[4] = -1e4 * y[2] - 6e7 * y[1], j[5] = -1e4 * y[1];
    j[6] = 0.0, j[7] = 6e7 * y[1], j[8] = 0.0;
    return 0;
}

static int refuses(double t, const double *y, double *dydt, void *context) {
    (void)t, (void)y, (void)context;
    dydt[0] = 0.0;
    return -1;
}

int main(void) {
    double y = 1.0, floor = 1e-12, step = 0.0;
    pl_rodas_options o = {.size = 1, .rhs = decay, .jacobian = decay_jacobian, .tolerance = 1e-9, .floor = &floor,
                          .smallest = 1e-12};
    pl_ode_statistics s = {0};
    CHECK(pl_rodas_advance(&o, 0.0, &y, 3.0, &step, &s, NULL) == PL_OK);
    CHECK_CLOSE(y, exp(-6.0), 1e-6);
    CHECK(s.steps > 0 && step > 0.0);

    /* Robertson's problem with Genesis's settings (src/core/rosenbrock.c at commit a7c127f, tolerance 1e-11, floors
     * 1e-14, 1e-18, 1e-14): the same algorithm, so the same numbers to rounding. */
    double r[3] = {1.0, 0.0, 0.0}, floors[3] = {1e-14, 1e-18, 1e-14}, t = 0.0;
    const double genesis[3][4] = {{0.4, 0.9851721138610029, 3.386395378970484e-05, 0.014794022185220277},
                                  {40.0, 0.71582706871732193, 9.1855347645162122e-06, 0.28416374574559122},
                                  {4e5, 0.004938274521044562, 1.9849940879853399e-08, 0.99506170562660234}};
    pl_rodas_options ro = {.size = 3, .rhs = robertson, .jacobian = robertson_jacobian, .tolerance = 1e-11,
                           .floor = floors, .smallest = 1e-30};
    step = 0.0;
    for (int k = 0; k < 3; k++) {
        CHECK(pl_rodas_advance(&ro, t, r, genesis[k][0] - t, &step, NULL, NULL) == PL_OK);
        t = genesis[k][0];
        for (int i = 0; i < 3; i++) CHECK_CLOSE(r[i], genesis[k][1 + i], 1e-12);
    }

    /* A callback that refuses: the model's failure, and y as it was. */
    double z = 2.0;
    pl_diagnostic d;
    pl_rodas_options bad = o;
    bad.rhs = refuses;
    step = 0.0;
    CHECK(pl_rodas_advance(&bad, 0.0, &z, 1.0, &step, NULL, &d) == PL_ERROR_MODEL && z == 2.0);
    CHECK(strstr(d.message, "callback") != NULL);
    pl_rodas_options short_of = o;
    short_of.max_trials = 1;
    z = 2.0, step = 1e-6;
    CHECK(pl_rodas_advance(&short_of, 0.0, &z, 10.0, &step, NULL, &d) == PL_ERROR_SOLVER && z == 2.0);
    return check_report("rodas");
}
