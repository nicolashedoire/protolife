/* Reaction networks (protolife/network.h): the text format, formulas, the balance and the rate law. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "protolife/protolife.h"

static void test_oregonator(void) {
    pl_network *net = NULL;
    pl_diagnostic d;
    CHECK(pl_network_load(PROTOLIFE_DATA "/networks/oregonator-reversible.net", &net, &d) == PL_OK);
    if (net == NULL) return;
    CHECK(pl_network_species_count(net) == 7);
    CHECK(pl_network_reaction_count(net) == 10);
    CHECK(pl_network_temperature(net) == 298.15);
    int x = pl_network_find(net, "X"), b = pl_network_find(net, "B"), h = pl_network_find(net, "Hplus");
    CHECK(x == 0 && pl_network_species_role(net, x) == PL_DYNAMIC);
    CHECK(pl_network_species_role(net, b) == PL_RESERVOIR && pl_network_species_molar(net, b) == 0.06);
    CHECK(pl_network_species_role(net, h) == PL_IMPLICIT && pl_network_species_charge(net, h) == 1);
    CHECK(strcmp(pl_network_species_formula(net, x), "HBrO2") == 0);
    CHECK(pl_network_reaction_constant(net, 2) == 2.40e6);
    CHECK(pl_network_reaction_lumped(net, 4) && !pl_network_reaction_lumped(net, 0));
    CHECK(strstr(pl_network_reaction_source(net, 0), "Dutt 2011") != NULL);
    CHECK(pl_network_find(net, "nothing") == -1);
    CHECK(pl_network_check_balance(net, &d) == PL_OK);
    CHECK(!pl_network_conserves(net, "Br"));  /* bromate is held from outside */
    CHECK(!pl_network_conserves(net, "Ce"));  /* the lumped steps do not keep cerium */
    pl_network_free(net);
    CHECK(pl_network_load(PROTOLIFE_DATA "/networks/e0-carbons.net", &net, &d) == PL_OK);
    CHECK(pl_network_conserves(net, "C") && !pl_network_conserves(net, "O") && !pl_network_conserves(net, "Xx"));
    pl_network_free(net);
}

static void test_balance(void) {
    pl_network *net = NULL;
    pl_diagnostic d;
    const char *text = "species A CH2O 0 dynamic 1\n"
                       "species B C2H4O2 0 dynamic 0\n"
                       "species Hplus H 1 implicit 1e-7\n"
                       "reaction 2 A -> B ; k=1 ; source=a dimerisation\n"
                       "reaction A -> B ; k=1 ; source=not balanced\n";
    CHECK(pl_network_parse(text, &net, &d) == PL_OK);
    CHECK(pl_network_check_balance(net, &d) == PL_ERROR_CONSERVATION);
    CHECK(d.code == 1 && strstr(d.message, "reaction 2") != NULL && strstr(d.message, "2 H on the left, 4 on the right"));
    pl_network_free(net);
    text = "species A H 0 dynamic 1\nspecies Hplus H 1 implicit 1\nreaction A -> ; k=1 ; implicit= -> Hplus ; source=s\n";
    CHECK(pl_network_parse(text, &net, &d) == PL_OK);
    CHECK(pl_network_check_balance(net, &d) == PL_ERROR_CONSERVATION && strstr(d.message, "charge 0 on the left"));
    pl_network_free(net);
    text = "species A CH2(COOH)2 0 dynamic 1\nspecies B C3H4O4 0 dynamic 0\nreaction A -> B ; k=1 ; source=s\n";
    CHECK(pl_network_parse(text, &net, &d) == PL_OK);
    CHECK(pl_network_check_balance(net, &d) == PL_OK);  /* malonic acid, written two ways */
    pl_network_free(net);
}

static void test_reading_errors(void) {
    pl_network *net = (pl_network *)&net;
    pl_diagnostic d;
    const char *wrong[][2] = {
        {"species A Xq 0 dynamic 1\n", "formula"},
        {"species A H 0 dynamic\n", "expected"},
        {"species A H 0 stuck 1\n", "role"},
        {"species A H 0 dynamic 1\nspecies A H 0 dynamic 1\n", "already exists"},
        {"species 2A H 0 dynamic 1\n", "digit"},
        {"species A H 0 dynamic 1\nreaction A -> C ; k=1 ; source=s\n", "unknown species C"},
        {"species A H 0 dynamic 1\nreaction A -> A ; source=s\n", "no k="},
        {"species A H 0 dynamic 1\nreaction A -> A ; k=1\n", "no source"},
        {"species A H 0 dynamic 1\nreaction A A ; k=1 ; source=s\n", "->"},
        {"species A H 0 dynamic 1\nreaction A -> A ; k=-1 ; source=s\n", ">= 0"},
        {"species A H 0 dynamic 1\nreaction A -> A ; k=one ; source=s\n", "not a number"},
        {"molecule A\n", "starts with"},
    };
    for (size_t i = 0; i < sizeof wrong / sizeof *wrong; i++) {
        pl_status status = pl_network_parse(wrong[i][0], &net, &d);
        CHECK(status == PL_ERROR_INPUT && net == NULL);
        if (strstr(d.message, wrong[i][1]) == NULL) fprintf(stderr, "case %zu: %s\n", i, d.message);
        CHECK(strstr(d.message, wrong[i][1]) != NULL);
    }
    CHECK(pl_network_parse("# a comment\n\nspecies A H 0 dynamic 1\nwhat\n", &net, &d) == PL_ERROR_INPUT);
    CHECK(d.code == 4 && strncmp(d.message, "line 4:", 7) == 0);
    CHECK(pl_network_load(PROTOLIFE_DATA "/networks/missing.net", &net, &d) == PL_ERROR_INPUT);
}

/* The Jacobian against central differences of the derivative, on the Oregonator at a state of its cycle. */
static void test_jacobian(void) {
    pl_network *net = NULL;
    CHECK(pl_network_load(PROTOLIFE_DATA "/networks/oregonator-reversible.net", &net, NULL) == PL_OK);
    if (net == NULL) return;
    int n = pl_network_species_count(net);
    double molar[7], up[7], down[7], dup[7], ddown[7], jacobian[49];
    for (int s = 0; s < n; s++) molar[s] = pl_network_species_molar(net, s);
    molar[0] = 4.6e-5, molar[1] = 4.6e-7, molar[2] = 5.7e-5;
    pl_network_jacobian(net, molar, jacobian);
    for (int q = 0; q < 3; q++) {
        double h = 1e-6 * molar[q];
        memcpy(up, molar, sizeof up), memcpy(down, molar, sizeof down);
        up[q] += h, down[q] -= h;
        pl_network_derivative(net, up, dup), pl_network_derivative(net, down, ddown);
        for (int s = 0; s < n; s++) {
            double estimate = (dup[s] - ddown[s]) / (2.0 * h);
            if (fabs(estimate) > 1e-12) CHECK_CLOSE(jacobian[s * n + q], estimate, 1e-6);
            else CHECK(fabs(jacobian[s * n + q]) < 1e-9);
        }
    }
    double rate[10], derivative[7];
    pl_network_rates(net, molar, rate);
    CHECK_CLOSE(rate[2], 2.40e6 * molar[0] * molar[1], 1e-15);  /* X + Y -> 2 P */
    pl_network_derivative(net, molar, derivative);
    CHECK(derivative[3] == 0.0 && derivative[5] == 0.0);        /* reservoir and implicit stay */
    pl_network_free(net);
}

static void test_building(void) {
    pl_network *net = NULL;
    pl_diagnostic d;
    CHECK(pl_network_create(&net, &d) == PL_OK);
    CHECK(pl_network_add_species(net, "glycolaldehyde", "C2H4O2", 0, PL_DYNAMIC, 0.01, &d) == PL_OK);
    CHECK(pl_network_add_species(net, "formaldehyde", "CH2O", 0, PL_RESERVOIR, 0.5, &d) == PL_OK);
    CHECK(pl_network_add_species(net, "glyceraldehyde", "C3H6O3", 0, PL_DYNAMIC, 0.0, &d) == PL_OK);
    CHECK(pl_network_add_reaction(net, "glycolaldehyde + formaldehyde -> glyceraldehyde", 1e-3, NULL, 0,
                                  "an aldol addition, for the test", &d) == PL_OK);
    CHECK(pl_network_add_reaction(net, "glycolaldehyde -> glyceraldehyde", 1.0, NULL, 0, "", &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_network_add_species(net, "bad name", "H", 0, PL_DYNAMIC, 0.0, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_network_add_species(net, "negative", "H", 0, PL_DYNAMIC, -1.0, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_network_set_temperature(net, 333.15, &d) == PL_OK && pl_network_temperature(net) == 333.15);
    CHECK(pl_network_check_balance(net, &d) == PL_OK);
    CHECK(pl_network_reaction_count(net) == 1 && pl_network_species_count(net) == 3);
    pl_network_free(net);
}

int main(void) {
    test_oregonator();
    test_balance();
    test_reading_errors();
    test_jacobian();
    test_building();
    return check_report("network");
}
