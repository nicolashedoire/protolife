/* The stiff integrator (protolife/ode.h): exact solutions, Genesis's reference, events, and honest failures. */
#include <math.h>
#include <string.h>

#include "check.h"
#include "protolife/protolife.h"

/* y' = -2 y, with its Jacobian: y = exp(-2 t). */
static int decay(double t, const double *y, double *dydt, void *context) {
    (void)t, (void)context;
    dydt[0] = -2.0 * y[0];
    return 0;
}

static int decay_jacobian(double t, const double *y, double *jacobian, void *context) {
    (void)t, (void)y, (void)context;
    jacobian[0] = -2.0;
    return 0;
}

static void test_decay(void) {
    double y0 = 1.0, absolute = 1e-12;
    pl_ode_options o = {.size = 1, .rhs = decay, .jacobian = decay_jacobian, .relative_tolerance = 1e-9,
                        .absolute_tolerance = &absolute};
    pl_ode *ode = NULL;
    pl_diagnostic d;
    CHECK(pl_ode_create(&o, 0.0, &y0, &ode, &d) == PL_OK);
    for (int k = 1; k <= 5; k++) {
        CHECK(pl_ode_advance(ode, k, NULL, &d) == PL_OK);
        CHECK(pl_ode_time(ode) == k);
        CHECK_CLOSE(pl_ode_state(ode)[0], exp(-2.0 * k), 1e-6);
    }
    pl_ode_free(ode);
}

/* A stiff scalar problem with an exact solution: u' = -1000 (u - cos t) - sin t, u(0) = 1, so u = cos t. */
static int stiff(double t, const double *y, double *dydt, void *context) {
    (void)context;
    dydt[0] = -1000.0 * (y[0] - cos(t)) - sin(t);
    return 0;
}

static void test_stiff(void) {
    double y0 = 1.0, absolute = 1e-12;
    pl_ode_options o = {.size = 1, .rhs = stiff, .relative_tolerance = 1e-9, .absolute_tolerance = &absolute};
    pl_ode *ode = NULL;
    CHECK(pl_ode_create(&o, 0.0, &y0, &ode, NULL) == PL_OK);
    CHECK(pl_ode_advance(ode, 10.0, NULL, NULL) == PL_OK);
    CHECK(fabs(pl_ode_state(ode)[0] - cos(10.0)) < 1e-7);
    pl_ode_statistics s;
    pl_ode_get_statistics(ode, &s);
    CHECK(s.steps > 0 && s.rhs_evaluations > s.steps);
    pl_ode_free(ode);
}

/* Robertson's problem (1966), the classic stiff chemistry: A -> B (0.04), 2 B -> B + C (3e7), B + C -> A + C (1e4). */
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

/* The same problem through Genesis's RODAS3 integrator (genesis src/core/rosenbrock.c at commit a7c127f, relative
 * tolerance 1e-11, floors 1e-14, 1e-18, 1e-14): time, then the three concentrations. */
static const double GENESIS_ROBERTSON[3][4] = {
    {0.4, 0.9851721138610029, 3.386395378970484e-05, 0.014794022185220277},
    {40.0, 0.71582706871732193, 9.1855347645162122e-06, 0.28416374574559122},
    {4e5, 0.004938274521044562, 1.9849940879853399e-08, 0.99506170562660234},
};

static void test_robertson(void) {
    double y0[3] = {1.0, 0.0, 0.0}, absolute[3] = {1e-14, 1e-18, 1e-14};
    pl_ode_options o = {.size = 3, .rhs = robertson, .jacobian = robertson_jacobian, .relative_tolerance = 1e-10,
                        .absolute_tolerance = absolute, .nonnegative = 1, .max_steps = 100000};
    pl_ode *ode = NULL;
    CHECK(pl_ode_create(&o, 0.0, y0, &ode, NULL) == PL_OK);
    for (int k = 0; k < 3; k++) {
        pl_diagnostic d;
        CHECK(pl_ode_advance(ode, GENESIS_ROBERTSON[k][0], NULL, &d) == PL_OK);
        const double *y = pl_ode_state(ode);
        for (int i = 0; i < 3; i++) CHECK_CLOSE(y[i], GENESIS_ROBERTSON[k][1 + i], 1e-7);  /* 2.4e-9 on macOS arm64 */
        CHECK(fabs(y[0] + y[1] + y[2] - 1.0) < 1e-12);
    }
    pl_ode_free(ode);
}

/* y' = 1 from 0, with the event y = 0.5. */
static int ramp(double t, const double *y, double *dydt, void *context) {
    (void)t, (void)y, (void)context;
    dydt[0] = 1.0;
    return 0;
}

static int half(double t, const double *y, double *g, void *context) {
    (void)t, (void)context;
    g[0] = y[0] - 0.5;
    return 0;
}

static void test_event(void) {
    double y0 = 0.0, absolute = 1e-12;
    pl_ode_options o = {.size = 1, .rhs = ramp, .relative_tolerance = 1e-10, .absolute_tolerance = &absolute,
                        .events = 1, .event = half};
    pl_ode *ode = NULL;
    int fired = 0;
    CHECK(pl_ode_create(&o, 0.0, &y0, &ode, NULL) == PL_OK);
    CHECK(pl_ode_advance(ode, 1.0, NULL, NULL) == PL_ERROR_ARGUMENT);
    CHECK(pl_ode_advance(ode, 1.0, &fired, NULL) == PL_OK);
    CHECK(fired == 1);
    CHECK(fabs(pl_ode_time(ode) - 0.5) < 1e-9);
    CHECK(fabs(pl_ode_state(ode)[0] - 0.5) < 1e-9);
    CHECK(pl_ode_advance(ode, 1.0, &fired, NULL) == PL_OK);
    CHECK(fired == 0 && pl_ode_time(ode) == 1.0);
    double y = 0.25;
    CHECK(pl_ode_restart(ode, 2.0, &y, NULL) == PL_OK);
    CHECK(pl_ode_advance(ode, 3.0, &fired, NULL) == PL_OK);
    CHECK(fired == 1 && fabs(pl_ode_time(ode) - 2.25) < 1e-9);
    pl_ode_free(ode);
}

/* y' = 1 until t = 0.5, then a right-hand side that is not finite. */
static int breaks(double t, const double *y, double *dydt, void *context) {
    (void)y, (void)context;
    dydt[0] = t <= 0.5 ? 1.0 : NAN;
    return 0;
}

static void test_model_failure(void) {
    double y0 = 0.0, absolute = 1e-12;
    pl_ode_options o = {.size = 1, .rhs = breaks, .relative_tolerance = 1e-8, .absolute_tolerance = &absolute};
    pl_ode *ode = NULL;
    pl_diagnostic d;
    CHECK(pl_ode_create(&o, 0.0, &y0, &ode, NULL) == PL_OK);
    CHECK(pl_ode_advance(ode, 1.0, NULL, &d) == PL_ERROR_MODEL);
    CHECK(d.status == PL_ERROR_MODEL && d.code < 0);
    CHECK(strstr(d.message, "not finite") != NULL);
    CHECK(pl_ode_time(ode) == 0.0 && pl_ode_state(ode)[0] == 0.0);  /* untouched */
    CHECK(pl_ode_advance(ode, 0.4, NULL, &d) == PL_OK);              /* and still usable */
    CHECK(fabs(pl_ode_state(ode)[0] - 0.4) < 1e-9);
    CHECK(d.status == PL_OK && d.message[0] == '\0');
    pl_ode_free(ode);
}

static void test_too_much_work(void) {
    double y0[3] = {1.0, 0.0, 0.0}, absolute[3] = {1e-14, 1e-18, 1e-14};
    pl_ode_options o = {.size = 3, .rhs = robertson, .relative_tolerance = 1e-10, .absolute_tolerance = absolute,
                        .max_steps = 5};
    pl_ode *ode = NULL;
    pl_diagnostic d;
    CHECK(pl_ode_create(&o, 0.0, y0, &ode, NULL) == PL_OK);
    CHECK(pl_ode_advance(ode, 4e5, NULL, &d) == PL_ERROR_SOLVER);
    CHECK(strstr(d.message, "CV_TOO_MUCH_WORK") != NULL);
    CHECK(d.time > 0.0 && d.time < 4e5);
    CHECK(pl_ode_time(ode) == 0.0 && pl_ode_state(ode)[0] == 1.0 && pl_ode_state(ode)[1] == 0.0);
    pl_ode_free(ode);
}

/* y' = -1 from 1 must cross zero at t = 1: with `nonnegative` the integration refuses, nothing is clipped. */
static int fall(double t, const double *y, double *dydt, void *context) {
    (void)t, (void)y, (void)context;
    dydt[0] = -1.0;
    return 0;
}

static void test_nonnegative(void) {
    double y0 = 1.0, absolute = 1e-12;
    pl_ode_options o = {.size = 1, .rhs = fall, .relative_tolerance = 1e-8, .absolute_tolerance = &absolute,
                        .nonnegative = 1};
    pl_ode *ode = NULL;
    pl_diagnostic d;
    CHECK(pl_ode_create(&o, 0.0, &y0, &ode, NULL) == PL_OK);
    CHECK(pl_ode_advance(ode, 0.5, NULL, NULL) == PL_OK);
    CHECK(fabs(pl_ode_state(ode)[0] - 0.5) < 1e-9);
    CHECK(pl_ode_advance(ode, 2.0, NULL, &d) == PL_ERROR_SOLVER);
    CHECK(pl_ode_time(ode) == 0.5 && fabs(pl_ode_state(ode)[0] - 0.5) < 1e-9);
    double negative = -1e-3;
    CHECK(pl_ode_restart(ode, 1.0, &negative, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_ode_time(ode) == 0.5);
    pl_ode_free(ode);
}

static void test_arguments(void) {
    double y0 = 1.0, absolute = 1e-12, zero = 0.0;
    pl_ode *ode = NULL;
    pl_diagnostic d;
    pl_ode_options o = {.size = 1, .rhs = decay, .relative_tolerance = 1e-8, .absolute_tolerance = &absolute};
    pl_ode_options empty = o, no_rhs = o, no_tolerance = o, bad_absolute = o, no_events = o;
    empty.size = 0, no_rhs.rhs = NULL, no_tolerance.relative_tolerance = 0.0, bad_absolute.absolute_tolerance = &zero;
    no_events.events = 1;
    const pl_ode_options *wrong[] = {&empty, &no_rhs, &no_tolerance, &bad_absolute, &no_events};
    for (size_t i = 0; i < sizeof wrong / sizeof *wrong; i++) {
        CHECK(pl_ode_create(wrong[i], 0.0, &y0, &ode, &d) == PL_ERROR_ARGUMENT);
        CHECK(ode == NULL && d.message[0] != '\0');
    }
    double nan_state = NAN;
    CHECK(pl_ode_create(&o, 0.0, &nan_state, &ode, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_ode_create(&o, 0.0, &y0, &ode, &d) == PL_OK);
    CHECK(pl_ode_advance(ode, 0.0, NULL, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_ode_advance(ode, -1.0, NULL, &d) == PL_ERROR_ARGUMENT);
    CHECK(strcmp(pl_status_name(d.status), "PL_ERROR_ARGUMENT") == 0);
    pl_ode_free(ode);
    pl_ode_free(NULL);
}

int main(void) {
    test_decay();
    test_stiff();
    test_robertson();
    test_event();
    test_model_failure();
    test_too_much_work();
    test_nonnegative();
    test_arguments();
    return check_report("ode");
}
