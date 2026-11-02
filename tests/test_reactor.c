/* The well-stirred reactor (protolife/reactor.h): exact solutions, an open tank, and Genesis's Oregonator. */
#include <math.h>
#include <string.h>

#include "check.h"
#include "protolife/protolife.h"

/* A <-> B, closed: A(t) = A_eq + (A0 - A_eq) exp(-(kf + kb) t), with A_eq = kb (A0 + B0) / (kf + kb). */
static void test_reversible(void) {
    pl_network *net = NULL;
    const char *text = "species A C4H8O4 0 dynamic 1e-2\n"
                       "species B C4H8O4 0 dynamic 0\n"
                       "reaction A -> B ; k=0.3 ; source=test\n"
                       "reaction B -> A ; k=0.1 ; source=test\n";
    CHECK(pl_network_parse(text, &net, NULL) == PL_OK);
    pl_reactor *r = NULL;
    pl_reactor_options o = {.relative_tolerance = 1e-10, .absolute_tolerance = 1e-16};
    CHECK(pl_reactor_create(net, &o, &r, NULL) == PL_OK);
    double equilibrium = 0.1 * 1e-2 / 0.4;
    for (int k = 1; k <= 20; k++) {
        CHECK(pl_reactor_advance(r, 0.5 * k, NULL) == PL_OK);
        const double *m = pl_reactor_molar(r);
        CHECK_CLOSE(m[0], equilibrium + (1e-2 - equilibrium) * exp(-0.4 * 0.5 * k), 1e-7);
        CHECK(fabs(m[0] + m[1] - 1e-2) < 1e-15);
    }
    pl_reactor_free(r);
    pl_network_free(net);
}

/* An open tank: A -> B at k, every species leaving at D, A fed at A_in. A(t) = A_s + (A0 - A_s) exp(-(D + k) t),
 * A_s = D A_in / (D + k). A reservoir catalyst C stays exactly at its value and scales k. */
static void test_open_tank(void) {
    pl_network *net = NULL;
    const char *text = "species A CH2O 0 dynamic 0\n"
                       "species B CH2O 0 dynamic 0\n"
                       "species C Fe 2 reservoir 2\n"
                       "reaction A + C -> B + C ; k=0.25 ; source=test\n";
    CHECK(pl_network_parse(text, &net, NULL) == PL_OK);
    double feed[3] = {1e-3, 0.0, 0.0};
    pl_reactor_options o = {.relative_tolerance = 1e-10, .absolute_tolerance = 1e-16, .dilution = 0.2, .feed = feed};
    pl_reactor *r = NULL;
    CHECK(pl_reactor_create(net, &o, &r, NULL) == PL_OK);
    double k = 0.25 * 2.0, steady = 0.2 * 1e-3 / (0.2 + k);
    for (int i = 1; i <= 10; i++) {
        CHECK(pl_reactor_advance(r, 2.0 * i, NULL) == PL_OK);
        CHECK_CLOSE(pl_reactor_molar(r)[0], steady * (1.0 - exp(-(0.2 + k) * 2.0 * i)), 1e-7);
        CHECK(pl_reactor_molar(r)[2] == 2.0);
    }
    CHECK(pl_reactor_set_molar(r, 2, 0.0, NULL) == PL_OK);  /* no catalyst: A only flows */
    CHECK(pl_reactor_advance(r, 200.0, NULL) == PL_OK);
    CHECK_CLOSE(pl_reactor_molar(r)[0], 1e-3, 1e-6);
    pl_reactor_free(r);
    pl_network_free(net);
}

/* The reversible Oregonator through Genesis's mean field (genesis tools: mean_field_advance, stiff, RODAS3, tolerance
 * 1e-11, floor 1e-16, commit a7c127f): time, then X, Y, Z in mol/L. It oscillates, so the later the time the more a
 * small phase difference shows. */
static const double GENESIS_OREGONATOR[4][4] = {
    {1.0, 4.6307103797956208e-07, 1.8849137118626269e-07, 2.136305535482423e-06},
    {5.0, 4.6715342042805072e-05, 4.645336463454819e-07, 5.6831581236006955e-05},
    {20.0, 3.8648178252789691e-08, 4.8747723514432317e-06, 4.9507033953436097e-05},
    {60.0, 4.5597414206831027e-08, 2.8077930942604967e-06, 2.8562087320766884e-05},
};

static void test_oregonator(void) {
    pl_network *net = NULL;
    CHECK(pl_network_load(PROTOLIFE_DATA "/networks/oregonator-reversible.net", &net, NULL) == PL_OK);
    if (net == NULL) return;
    pl_reactor_options o = {.relative_tolerance = 1e-11, .absolute_tolerance = 1e-18, .max_steps = 1000000};
    pl_reactor *r = NULL;
    pl_diagnostic d;
    CHECK(pl_reactor_create(net, &o, &r, &d) == PL_OK);
    for (int k = 0; k < 4; k++) {
        CHECK(pl_reactor_advance(r, GENESIS_OREGONATOR[k][0], &d) == PL_OK);
        for (int s = 0; s < 3; s++)  /* 6.2e-10 at most on macOS arm64 */
            CHECK_CLOSE(pl_reactor_molar(r)[s], GENESIS_OREGONATOR[k][1 + s], 1e-7);
    }
    CHECK(pl_reactor_molar(r)[3] == 0.06 && pl_reactor_molar(r)[4] == 1e-4);
    pl_reactor_free(r);
    pl_network_free(net);
}

static void test_arguments(void) {
    pl_network *net = NULL;
    pl_reactor *r = NULL;
    pl_diagnostic d;
    CHECK(pl_network_parse("species W H2O 0 reservoir 55.5\n", &net, NULL) == PL_OK);
    CHECK(pl_reactor_create(net, NULL, &r, &d) == PL_ERROR_ARGUMENT && r == NULL);
    pl_network_free(net);
    CHECK(pl_network_parse("species A H 0 dynamic 1\nspecies H H 1 implicit 1\n", &net, NULL) == PL_OK);
    pl_reactor_options negative = {.dilution = -1.0};
    CHECK(pl_reactor_create(net, &negative, &r, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_reactor_create(net, NULL, &r, &d) == PL_OK);
    CHECK(pl_reactor_set_molar(r, 1, 2.0, &d) == PL_ERROR_ARGUMENT);   /* implicit */
    CHECK(pl_reactor_set_molar(r, 0, -1.0, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_reactor_advance(r, 1.0, &d) == PL_OK && pl_reactor_molar(r)[0] == 1.0);
    pl_reactor_free(r);
    pl_network_free(net);
}

int main(void) {
    test_reversible();
    test_open_tank();
    test_oregonator();
    test_arguments();
    return check_report("reactor");
}
