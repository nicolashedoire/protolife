/* Fatty acid vesicles (protolife/vesicles.h): the checks of Genesis's tests/test_vesicles.c, through the library. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "protolife/protolife.h"

#define PI 3.14159265358979323846
#define AVOGADRO 6.02214076e23
#define GAS 8.314462618

/* A medium of 0.3 osm/L of a kept solute, and one record of 100 nm oleate vesicles holding 80 % of their sphere's
 * worth of osmoles, plus `extra_osm` osm/m^3. */
static pl_vesicles *population(double extra_osm, int *kept) {
    pl_vesicle_constants k = pl_vesicle_oleate();
    pl_vesicles *p = NULL;
    CHECK(pl_vesicles_create(&k, 0, &p, NULL) == PL_OK);
    CHECK(pl_vesicles_add_solute(p, "buffer", 0.0, 1.0, kept, NULL) == PL_OK);
    double radius = 50e-9, molecules = pl_vesicle_molecules_for_radius(&k, radius);
    double sphere = 4.0 / 3.0 * PI * radius * radius * radius, inside = (300.0 + extra_osm) * 0.8 * sphere;
    CHECK(pl_vesicles_add(p, 3e-3 * AVOGADRO / molecules, molecules, &inside, 7, NULL, NULL) == PL_OK);
    CHECK(pl_vesicles_set_outside(p, *kept, 0.3, NULL) == PL_OK);
    CHECK(pl_vesicles_set_monomer(p, k.cvc * pl_vesicles_medium(p, NULL), NULL) == PL_OK);
    return p;
}

static void test_shapes(void) {
    pl_vesicle_constants k = pl_vesicle_oleate();
    pl_vesicle_shape loose = pl_vesicle_shape_of(&k, pl_vesicle_molecules_for_radius(&k, 50e-9), 1e-19, 300.0);
    CHECK_CLOSE(loose.volume, 1e-19 / 300.0, 1e-12);
    CHECK(loose.tension == 0.0 && !loose.spherical);
    double r = 50e-9, sphere = 4.0 / 3.0 * PI * r * r * r;
    pl_vesicle_shape round = pl_vesicle_shape_of(&k, pl_vesicle_molecules_for_radius(&k, r), 400.0 * sphere, 300.0);
    CHECK_CLOSE(round.tension, 100.0 * GAS * k.kelvin * r / 2.0, 1e-12);  /* Laplace */
}

static void test_balance(void) {
    int kept;
    pl_vesicles *p = population(0.0, &kept);
    double acid = pl_vesicles_fatty_acid(p), solute = pl_vesicles_solute_total(p, kept);
    pl_vesicle_shape before = pl_vesicles_shape(p, 0);
    CHECK(pl_vesicles_advance(p, 60.0, NULL) == PL_OK);
    CHECK_CLOSE(pl_vesicles_shape(p, 0).area, before.area, 1e-9);  /* in balance with its monomer, it stays */
    CHECK_CLOSE(pl_vesicles_fatty_acid(p), acid, 1e-12);
    CHECK_CLOSE(pl_vesicles_solute_total(p, kept), solute, 1e-12);
    CHECK(pl_vesicles_time(p) == 60.0 && pl_vesicles_find(p, 7) == 0);
    pl_vesicles_free(p);
}

static void test_rupture(void) {
    int kept;
    pl_vesicles *p = population(2000.0, &kept);  /* 2 osm/L too many: far past rupture */
    double solute = pl_vesicles_solute_total(p, kept);
    CHECK(pl_vesicles_advance(p, 1.0, NULL) == PL_OK);
    pl_vesicle_ledger ledger = pl_vesicles_ledger(p);
    CHECK(ledger.ruptures == 1 && ledger.leaked[kept] > 0.0);
    CHECK(pl_vesicles_shape(p, 0).tension <= pl_vesicle_oleate().rupture * (1.0 + 1e-9));
    CHECK_CLOSE(pl_vesicles_solute_total(p, kept), solute, 1e-12);  /* what it lets out stays in the vessel */
    pl_vesicles_free(p);
}

static void test_growth_and_division(void) {
    int kept;
    pl_vesicles *p = population(0.0, &kept);
    pl_vesicle_shape before = pl_vesicles_shape(p, 0);
    double fed = pl_vesicles_membranes(p), acid = pl_vesicles_fatty_acid(p), solute = pl_vesicles_solute_total(p, kept);
    CHECK(pl_vesicles_set_shear(p, 20.0, NULL) == PL_OK);  /* flaccid unfed vesicles do not divide */
    CHECK(pl_vesicles_advance(p, 1.0, NULL) == PL_OK && pl_vesicles_ledger(p).divided == 0.0);
    CHECK(pl_vesicles_set_shear(p, 0.0, NULL) == PL_OK);
    pl_vesicles_feed_micelles(p, fed);
    CHECK_CLOSE(pl_vesicles_entering(p), fed, 1e-12);
    int grew = 1;
    for (int i = 0; i < 120 && grew; i++) grew = pl_vesicles_advance(p, 1.0, NULL) == PL_OK;
    CHECK(grew && pl_vesicles_records(p) == 2 && pl_vesicles_record(p, 1).origin == 1);
    CHECK(pl_vesicles_ledger(p).born > 0.0 && pl_vesicles_shape(p, 0).area > 1.7 * before.area);
    CHECK_CLOSE(pl_vesicles_fatty_acid(p), acid + fed, 1e-10);
    CHECK_CLOSE(pl_vesicles_solute_total(p, kept), solute, 1e-10);  /* what new vesicles enclose leaves the medium */
    pl_vesicle v = pl_vesicles_record(p, 0);
    double count = v.count, own = v.count * v.molecules;
    CHECK(pl_vesicles_set_shear(p, 20.0, NULL) == PL_OK);
    CHECK(pl_vesicles_advance(p, 1.0, NULL) == PL_OK);
    v = pl_vesicles_record(p, 0);
    pl_vesicle_shape daughter = pl_vesicles_shape(p, 0);
    double ratio = v.count / count;
    CHECK(ratio >= 2.0 && fabs(ratio - floor(ratio + 0.5)) < 1e-12 && v.generation == 1);
    CHECK(daughter.volume / (4.0 / 3.0 * PI * pow(daughter.radius, 3.0)) >= pl_vesicle_oleate().division_reduced);
    CHECK_CLOSE(v.count * v.molecules, own, 1e-6);  /* sharing their membrane */
    pl_vesicles_free(p);
}

/* Permeability by carbons, as Genesis's retention world: one or none fast, two to four `small`, five and six or more
 * from Sacerdote and Szostak 2005. */
static double by_carbons(const pl_network *network, int species, void *context) {
    int carbons = pl_network_species_atoms(network, species, "C");
    double small = *(const double *)context;
    return carbons <= 1 ? 1e-6 : carbons <= 4 ? small : carbons == 5 ? 1e-10 : 5e-13;
}

static double carbon(const pl_vesicles *p, const pl_network *net, double *renewed) {
    double total = 0.0, brought = 0.0;
    pl_vesicle_ledger ledger = pl_vesicles_ledger(p);
    for (int s = 0; s < pl_vesicles_solutes(p); s++) {
        int k = pl_vesicles_species_of(p, s), c = k >= 0 ? pl_network_species_atoms(net, k, "C") : 0;
        total += c * pl_vesicles_solute_total(p, s);
        brought += c * ledger.renewed[s];
    }
    *renewed = brought;
    return total;
}

/* An E0 core of Origin inside and outside, the medium renewed: every carbon is accounted for. */
static void test_chemistry_inside(void) {
    pl_network *net = NULL;
    pl_diagnostic d;
    CHECK(pl_network_load(PROTOLIFE_DATA "/networks/e0-random1.net", &net, &d) == PL_OK);
    if (net == NULL) return;
    pl_vesicle_constants k = pl_vesicle_oleate();
    k.kelvin = pl_network_temperature(net);
    pl_vesicles *p = NULL;
    double small = 2.9e-10;
    int buffer, n = pl_network_species_count(net), food = pl_vesicles_solute_of(p, 0);
    CHECK(pl_vesicles_create(&k, 0, &p, NULL) == PL_OK);
    CHECK(pl_vesicles_add_solute(p, "bicine", 0.0, 1.635, &buffer, NULL) == PL_OK);
    CHECK(pl_vesicles_attach_network(p, net, by_carbons, &small, &d) == PL_OK);
    food = pl_vesicles_solute_of(p, pl_network_find(net, "C=O"));
    CHECK(food >= 0 && pl_vesicles_solutes(p) == 1 + n - 1);  /* every species but water */
    double radius = 500e-9, molecules = pl_vesicle_molecules_for_radius(&k, radius);
    double volume = 4.0 / 3.0 * PI * pow(radius * sqrt(1.0 - 0.27), 3.0);
    double *seeded = calloc((size_t)pl_vesicles_solutes(p), sizeof *seeded);
    double *fresh = calloc((size_t)pl_vesicles_solutes(p), sizeof *fresh);
    seeded[buffer] = fresh[buffer] = 0.2 * 1000.0 * volume;
    for (int s = 0; s < pl_vesicles_solutes(p); s++) {
        int species = pl_vesicles_species_of(p, s);
        if (species >= 0) seeded[s] = pl_network_species_molar(net, species) * 1000.0 * volume;
    }
    fresh[food] = 0.5 * 1000.0 * volume;
    double count = 1.5e-3 * AVOGADRO / molecules;
    CHECK(pl_vesicles_add(p, count, molecules, seeded, 1, NULL, NULL) == PL_OK);
    CHECK(pl_vesicles_add(p, count, molecules, fresh, 2, NULL, NULL) == PL_OK);
    CHECK(pl_vesicles_add_solute(p, "late", 0.0, 1.0, NULL, NULL) == PL_ERROR_ARGUMENT);  /* solutes come first */
    CHECK(pl_vesicles_set_feed(p, buffer, 0.2, NULL) == PL_OK && pl_vesicles_set_feed(p, food, 0.5, NULL) == PL_OK);
    CHECK(pl_vesicles_set_monomer_feed(p, k.cvc, NULL) == PL_OK && pl_vesicles_set_renewal(p, 5e-2, 0, NULL) == PL_OK);
    for (int pass = 0; pass < 3; pass++) {  /* the medium is what the vesicles leave of the vessel */
        double medium = pl_vesicles_medium(p, NULL);
        pl_vesicles_set_outside(p, buffer, 0.2 * medium, NULL), pl_vesicles_set_outside(p, food, 0.5 * medium, NULL);
        pl_vesicles_set_monomer(p, k.cvc * medium, NULL);
    }
    double renewed, carbon0 = carbon(p, net, &renewed);
    int ran = 1;
    for (int i = 0; i < 6 && ran; i++) ran = pl_vesicles_advance(p, 10.0, &d) == PL_OK;
    if (!ran) fprintf(stderr, "%s\n", d.message);
    double after = carbon(p, net, &renewed);
    CHECK(ran);
    CHECK_CLOSE(after - renewed, carbon0, 1e-9);  /* carbon in, out and reacted closes */
    double sugar = 0.0;
    pl_vesicle fresh_record = pl_vesicles_record(p, pl_vesicles_find(p, 2));
    for (int s = 0; s < pl_vesicles_solutes(p); s++) {
        int species = pl_vesicles_species_of(p, s);
        if (species >= 0 && pl_network_species_atoms(net, species, "C") >= 2) sugar += fresh_record.inside[s];
    }
    CHECK(sugar > 0.0);  /* the fresh vesicles make sugars of their own */
    pl_vesicles_dilute(p, 0.5);  /* half the vessel replaced: what the vesicles held leaves with them */
    after = carbon(p, net, &renewed);
    CHECK_CLOSE(after - renewed, carbon0, 1e-9);
    /* the balances checked at every step: the fatty acid, the buffer, carbon; not hydrogen or oxygen (water) */
    const pl_vesicle_balance *b;
    int balances = pl_vesicles_balances(p, &b), found = 0;
    for (int i = 0; i < balances; i++) {
        found |= (strcmp(b[i].name, "fatty acid") == 0) | (strcmp(b[i].name, "bicine") == 0) << 1 |
                 (strcmp(b[i].name, "C") == 0) << 2;
        CHECK(strcmp(b[i].name, "H") != 0 && strcmp(b[i].name, "O") != 0);
        CHECK(fabs(b[i].drift) <= 1e-12 * fabs(b[i].value) + 1e-18);
    }
    CHECK(balances == 3 && found == 7);
    free(seeded), free(fresh);
    pl_vesicles_free(p);
    pl_network_free(net);
}

/* A limit no step can meet: the step is undone, the population as it was, and the balance named. */
static void test_conservation_refused(void) {
    int kept;
    pl_vesicles *p = population(0.0, &kept);
    pl_vesicles_feed_micelles(p, 0.5 * pl_vesicles_membranes(p));
    CHECK(pl_vesicles_advance(p, 1.0, NULL) == PL_OK);
    double acid = pl_vesicles_fatty_acid(p), time = pl_vesicles_time(p), monomer = pl_vesicles_monomer(p);
    pl_vesicle v = pl_vesicles_record(p, 0);
    int records = pl_vesicles_records(p);
    pl_diagnostic d;
    CHECK(pl_vesicles_set_conservation(p, 1e-300, &d) == PL_OK);
    CHECK(pl_vesicles_advance(p, 5.0, &d) == PL_ERROR_CONSERVATION);
    CHECK(strstr(d.message, "the balance of") != NULL && strstr(d.message, "undone") != NULL);
    CHECK(pl_vesicles_time(p) == time && pl_vesicles_fatty_acid(p) == acid && pl_vesicles_monomer(p) == monomer);
    CHECK(pl_vesicles_record(p, 0).molecules == v.molecules && pl_vesicles_records(p) == records);
    CHECK(pl_vesicles_set_conservation(p, 1e-9, &d) == PL_OK);
    CHECK(pl_vesicles_advance(p, 5.0, &d) == PL_OK);  /* and it goes on with a sound limit */
    const pl_vesicle_balance *b;
    CHECK(pl_vesicles_balances(p, &b) == 2 && strcmp(b[0].name, "fatty acid") == 0 && strcmp(b[1].name, "buffer") == 0);
    CHECK(pl_vesicles_ledger(p).fed > 0.0);
    pl_vesicles_free(p);
}

static void test_arguments(void) {
    pl_vesicle_constants k = pl_vesicle_oleate(), bad = k;
    bad.cvc = -1.0;
    pl_vesicles *p = NULL;
    pl_diagnostic d;
    CHECK(pl_vesicles_create(&bad, 0, &p, &d) == PL_ERROR_ARGUMENT && p == NULL);
    CHECK(pl_vesicles_create(&k, 1, &p, &d) == PL_OK);
    CHECK(pl_vesicles_add_solute(p, "salt", -1.0, 2.0, NULL, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_vesicles_add(p, 1e15, 1e5, NULL, 0, NULL, &d) == PL_OK);
    CHECK(pl_vesicles_add(p, 1e15, 1e5, NULL, 0, NULL, &d) == PL_ERROR_ARGUMENT);  /* one record at most */
    CHECK(pl_vesicles_advance(p, 0.0, &d) == PL_ERROR_ARGUMENT);
    CHECK(pl_vesicles_set_shear(p, -1.0, &d) == PL_ERROR_ARGUMENT);
    pl_vesicles_free(p);
}

int main(void) {
    test_shapes();
    test_balance();
    test_rupture();
    test_growth_and_division();
    test_chemistry_inside();
    test_conservation_refused();
    test_arguments();
    return check_report("vesicles");
}
