/* Taneja and Higgs's protocell (protolife/capsule.h): the results of their paper (Life 15, 2025), as regression tests.
 * Their standard set throughout; times in the units of their rate constants; step 0.01 as in their program. */
#include <math.h>
#include <string.h>

#include "check.h"
#include "protolife/protolife.h"

static pl_network *net;

typedef struct {
    double elip, rmem, rflip, muw;
    int equal, inactive;
} setting;

static pl_capsule *cell(setting s) {
    pl_capsule_parameters k = pl_capsule_taneja_higgs();
    k.elip = s.elip * k.cstar, k.rmem = s.rmem, k.rflip = s.rflip, k.equal_division = s.equal;
    double mu[3] = {0.01, 0.0, s.muw}, outside[3] = {1.0, 0.0, 0.0};  /* food in, catalyst kept, waste out */
    pl_capsule *c = NULL;
    CHECK(pl_capsule_create(&k, net, mu, outside, 5505, &c, NULL) == PL_OK);
    if (s.inactive) pl_capsule_set_molar(c, 0, 1.0, NULL), pl_capsule_set_molar(c, 1, 0.0, NULL);  /* no catalyst */
    CHECK(pl_capsule_settle(c, 1000.0, 0.01, NULL) == PL_OK);
    return c;
}

typedef struct {
    int divisions, bursts;
    double reduced, small;  /* means at division: v, and the smaller daughter's share */
} tally;

static tally count(const pl_capsule *c) {
    const pl_capsule_event *e;
    int n = pl_capsule_events(c, &e);
    tally t = {0};
    for (int i = 0; i < n; i++) {
        if (e[i].kind == PL_CAPSULE_BURST) {
            t.bursts++;
            continue;
        }
        t.divisions++, t.reduced += e[i].reduced, t.small += fmin(e[i].share, 1.0 - e[i].share);
    }
    if (t.divisions > 0) t.reduced /= t.divisions, t.small /= t.divisions;
    return t;
}

static const setting STANDARD = {1.01, 0.1, 0.1, 0.015, 0, 0};

/* The fixed compartment: the active steady state solves C1^2 (1 - C1) = 2 w^2 / (k mu1) = 0.02. */
static void test_chemistry(void) {
    pl_capsule *c = cell(STANDARD);
    double c1 = pl_capsule_molar(c, 0);
    CHECK_CLOSE(c1 * c1 * (1.0 - c1), 0.02, 1e-4);
    CHECK_CLOSE(pl_capsule_molar(c, 1), 0.01 * (1.0 - c1) / 0.02, 1e-4);  /* C2 = mu1 (1 - C1) / (2 w) */
    pl_capsule_free(c);
}

/* Homeostasis (their Fig. 3a): with Elip = C* the cell stops at s = 2 w (muW - mu1) / (mu1 muW) = 2/3, that is
 * A = 9/4 A0 and V = 27/8 V0. */
static void test_homeostasis(void) {
    setting s = STANDARD;
    s.elip = 1.0;
    pl_capsule *c = cell(s);
    CHECK(pl_capsule_advance(c, 20000.0, 0.01, NULL) == PL_OK);
    pl_capsule_shape shape = pl_capsule_get_shape(c);
    CHECK_CLOSE(shape.volume, 27.0 / 8.0, 1e-4);
    CHECK_CLOSE(shape.area, 9.0 / 4.0, 1e-4);
    CHECK(count(c).divisions == 0);  /* division needs Elip > C* */
    pl_capsule_free(c);
}

/* With muW = mu1 there is no stationary state (Fig. 3c): the cell keeps growing and the reaction collapses. */
static void test_no_stationary_state(void) {
    setting s = STANDARD;
    s.elip = 1.0, s.muw = 0.01;
    pl_capsule *c = cell(s);
    double catalyst = pl_capsule_molar(c, 1);
    CHECK(pl_capsule_advance(c, 20000.0, 0.01, NULL) == PL_OK);
    CHECK(pl_capsule_get_shape(c).volume > 100.0 && pl_capsule_molar(c, 1) < 0.1 * catalyst);
    pl_capsule_free(c);
}

/* The equal rule (Fig. 4): every division at v = 1/2, into halves. */
static void test_equal_division(void) {
    setting s = STANDARD;
    s.equal = 1;
    pl_capsule *c = cell(s);
    CHECK(pl_capsule_advance(c, 20000.0, 0.01, NULL) == PL_OK);
    tally t = count(c);
    CHECK(t.divisions >= 5 && fabs(t.reduced - 0.5) < 1e-3 && t.small == 0.5 && t.bursts == 0);
    pl_capsule_free(c);
}

/* The energy rule and the flip-flop (Fig. 6): divisions up to rflip = 0.18 and none past it; the slower the
 * flipping, the closer v is to 1 at division and the more unequal the daughters. */
static void test_flip_flop(void) {
    const double rates[] = {0.01, 0.1, 0.18, 0.19};
    tally t[4];
    for (int i = 0; i < 4; i++) {
        setting s = STANDARD;
        s.rflip = rates[i];
        pl_capsule *c = cell(s);
        CHECK(pl_capsule_advance(c, 20000.0, 0.01, NULL) == PL_OK);
        t[i] = count(c);
        pl_capsule_free(c);
    }
    CHECK(t[0].divisions > 0 && t[1].divisions > 0 && t[2].divisions > 0 && t[3].divisions == 0);
    CHECK(t[0].reduced > 0.9 && t[0].reduced > t[1].reduced && t[1].reduced > t[2].reduced && t[2].reduced > 0.5);
    CHECK(t[0].small < 0.05 && t[0].small < t[1].small);
}

/* Without lipid exchange the cell only bursts (Fig. 8a); with slow exchange it bursts, then divides (Fig. 8b). */
static void test_bursts(void) {
    setting s = STANDARD;
    s.rmem = s.rflip = 0.0;
    pl_capsule *c = cell(s);
    CHECK(pl_capsule_advance(c, 20000.0, 0.01, NULL) == PL_OK);
    tally t = count(c);
    CHECK(t.divisions == 0 && t.bursts > 50);
    pl_capsule_free(c);
    s.rmem = s.rflip = 0.01;
    c = cell(s);
    CHECK(pl_capsule_advance(c, 40000.0, 0.01, NULL) == PL_OK);
    t = count(c);
    CHECK(t.bursts > 0 && t.divisions > 0);
    pl_capsule_free(c);
}

/* Inactive cells (Fig. 7): at Elip = 1.01 C* the area grows at a nearly fixed volume and the cell never divides; at
 * 1.02 C* it divides a few times, its volume shrinking, then grows in area without dividing. Both end as thin tubes,
 * where the authors stop. */
static void test_inactive(void) {
    setting s = STANDARD;
    s.inactive = 1;
    pl_capsule *c = cell(s);
    pl_diagnostic d;
    CHECK(pl_capsule_advance(c, 20000.0, 0.01, &d) == PL_ERROR_MODEL && strstr(d.message, "tube") != NULL);
    pl_capsule_shape shape = pl_capsule_get_shape(c);
    CHECK(count(c).divisions == 0 && shape.lipid_area > 20.0 && shape.volume < 1.0);
    pl_capsule_free(c);
    s.elip = 1.02;
    c = cell(s);
    CHECK(pl_capsule_advance(c, 20000.0, 0.01, &d) == PL_ERROR_MODEL);
    tally t = count(c);
    shape = pl_capsule_get_shape(c);
    CHECK(t.divisions >= 1 && t.divisions <= 5 && shape.volume < 0.5 && shape.lipid_area > 2.0);
    pl_capsule_free(c);
}

int main(void) {
    CHECK(pl_network_load(PROTOLIFE_DATA "/networks/taneja-higgs.net", &net, NULL) == PL_OK);
    if (net == NULL) return check_report("capsule");
    test_chemistry();
    test_homeostasis();
    test_no_stationary_state();
    test_equal_division();
    test_flip_flop();
    test_bursts();
    test_inactive();
    pl_network_free(net);
    return check_report("capsule");
}
