/* The world of the reproduction example (reproduction.h): its defaults, its building and installation, its steps. The
 * tests drive it too (tests/test_selection.c). */
#include <math.h>
#include <stdlib.h>

#include "reproduction.h"

options defaults(void) {
    return (options){"data/networks/e0-random2.net", "recordings/step-6-random2", "random2", 1, 500e-9, 3e-3, 0.9,
                     2.9e-10, 0.5, 5e-2, 361, 10.0, 0.5, {4, 3600.0, 900.0, 0.5, 1.0, 1, 60.0, 20.0, 1}, 0, 0, 0.0, 1};
}

/* Permeability by carbons: one or none fast (assumed); two to four `small`, not measured through oleate; five
 * 1e-8 cm/s (aldopentoses, Sacerdote and Szostak 2005); six or more 5e-11 cm/s (glucose, same source). */
static double by_carbons(const pl_network *network, int species, void *context) {
    int carbons = pl_network_species_atoms(network, species, "C");
    double small = *(const double *)context;
    return carbons <= 1 ? 1e-6 : carbons <= 4 ? small : carbons == 5 ? 1e-10 : 5e-13;
}

/* Amounts (mol) in one vesicle of `volume` (m^3) holding the buffer and the network at `molar` (per species). */
static void filled(const world *w, const double *molar, double volume, double *inside) {
    for (int s = 0; s < pl_vesicles_solutes(w->p); s++) {
        int k = pl_vesicles_species_of(w->p, s);
        inside[s] = s == w->buffer ? BUFFER_MOLAR * 1000.0 * volume : k >= 0 ? molar[k] * 1000.0 * volume : 0.0;
    }
}

static pl_status seed(world *w, pl_diagnostic *d) {
    const options *o = &w->o;
    pl_vesicle_constants k = pl_vesicle_oleate();
    int n = pl_network_species_count(w->network), solutes = pl_vesicles_solutes(w->p);
    double *core = calloc((size_t)n, sizeof *core), *fresh = calloc((size_t)n, sizeof *fresh);
    double *inside = calloc((size_t)solutes, sizeof *inside);
    if (core == NULL || fresh == NULL || inside == NULL) return free(core), free(fresh), free(inside), PL_ERROR_MEMORY;
    for (int s = 0; s < n; s++) core[s] = pl_network_species_role(w->network, s) == PL_DYNAMIC ? pl_network_species_molar(w->network, s) : 0.0;
    fresh[pl_vesicles_species_of(w->p, w->food)] = o->food;
    double volume = 4.0 / 3.0 * PI * pow(o->radius * sqrt(1.0 - EXCESS), 3.0);
    double molecules = pl_vesicle_molecules_for_radius(&k, o->radius);
    filled(w, o->twins ? fresh : core, volume, inside);  /* twins: the seeded lineage holds fresh medium too */
    pl_status status = pl_vesicles_add(w->p, o->acid * o->seeded * AVOGADRO / molecules, molecules, inside, SEEDED, NULL, d);
    filled(w, fresh, volume, inside);
    if (status == PL_OK)
        status = pl_vesicles_add(w->p, o->acid * (1.0 - o->seeded) * AVOGADRO / molecules, molecules, inside, FRESH, NULL, d);
    free(core), free(fresh), free(inside);
    return status;
}

pl_status build(world *w, pl_diagnostic *d) {
    pl_status status = pl_network_load(w->o.network, &w->network, d);
    if (status != PL_OK) return status;
    pl_vesicle_constants k = pl_vesicle_oleate();
    k.kelvin = pl_network_temperature(w->network);
    if ((status = pl_vesicles_create(&k, 0, &w->p, d)) != PL_OK) return status;
    if ((status = pl_vesicles_add_solute(w->p, "bicine", 0.0, BUFFER_OSMOTIC, &w->buffer, d)) != PL_OK) return status;
    if ((status = pl_vesicles_attach_network(w->p, w->network, by_carbons, &w->o.small, d)) != PL_OK) return status;
    w->food = pl_vesicles_solute_of(w->p, pl_network_find(w->network, "C=O"));
    if (w->food < 0) return PL_ERROR_INPUT;
    if ((status = seed(w, d)) != PL_OK) return status;
    pl_vesicles_set_lifeless(w->p, !w->o.living);
    if (w->o.vesicles > 0.0) {  /* a bag small enough to hold that many seeded vesicles: its dilutions are drawn */
        double seeded = pl_vesicles_record(w->p, pl_vesicles_find(w->p, SEEDED)).count;
        if ((status = pl_vesicles_set_bag(w->p, w->o.vesicles / seeded, w->o.draw, d)) != PL_OK) return status;
    }
    pl_vesicles_set_feed(w->p, w->buffer, BUFFER_MOLAR, d), pl_vesicles_set_feed(w->p, w->food, w->o.food, d);
    pl_vesicles_set_monomer_feed(w->p, k.cvc, d);
    pl_vesicles_set_renewal(w->p, w->o.renewal, 1, d);  /* a dialysis bag: the fatty acid is held back */
    for (int pass = 0; pass < 3; pass++) {  /* the medium is what the vesicles leave of the vessel */
        double medium = pl_vesicles_medium(w->p, NULL);
        pl_vesicles_set_outside(w->p, w->buffer, BUFFER_MOLAR * medium, d);
        pl_vesicles_set_outside(w->p, w->food, w->o.food * medium, d);
        pl_vesicles_set_monomer(w->p, k.cvc * medium, d);
    }
    return PL_OK;
}

void dismantle(world *w) {
    pl_vesicles_free(w->p);
    pl_network_free(w->network);
    w->p = NULL, w->network = NULL;
}

/* Carbon in sugars (two carbons or more) held by one vesicle of a record (mol), or in the medium (mol/L) for -1. */
double sugar_carbon(const world *w, int record) {
    double carbon = 0.0;
    pl_vesicle v = record >= 0 ? pl_vesicles_record(w->p, record) : (pl_vesicle){0};
    for (int s = 0; s < pl_vesicles_solutes(w->p); s++) {
        int k = pl_vesicles_species_of(w->p, s), c = k >= 0 ? pl_network_species_atoms(w->network, k, "C") : 0;
        if (c >= 2) carbon += c * (record >= 0 ? v.inside[s] : pl_vesicles_outside(w->p, s));
    }
    return record >= 0 ? carbon : carbon / pl_vesicles_medium(w->p, NULL);
}

double fatty_acid_of(const world *w, int label) {
    double acid = 0.0;
    for (int r = 0; r < pl_vesicles_records(w->p); r++) {
        pl_vesicle v = pl_vesicles_record(w->p, r);
        if (v.label == label) acid += v.count * v.molecules / AVOGADRO;
    }
    return acid;
}

pl_status step(world *w, pl_diagnostic *d) {
    double before = pl_vesicles_time(w->p);
    pl_status status = pl_vesicles_advance(w->p, w->o.dt, d);
    if (status == PL_OK) pl_cycles_apply(w->p, &w->o.cycles, before, pl_vesicles_time(w->p));
    return status;
}

/* The installation, unrecorded, up to a step before the first cycle; there, the reference; then that step. */
pl_status install(world *w, pl_diagnostic *d) {
    pl_status status = PL_OK;
    while (status == PL_OK && pl_vesicles_time(w->p) + w->o.dt < w->o.cycles.start - w->o.dt / 2.0) status = step(w, d);
    int r = pl_vesicles_find(w->p, SEEDED);
    if (status == PL_OK && r < 0) return PL_ERROR_MODEL;  /* the seeded vesicles came apart */
    if (status != PL_OK) return status;
    w->count0 = pl_vesicles_record(w->p, r).count, w->held0 = sugar_carbon(w, r);
    pl_vesicle_ledger ledger = pl_vesicles_ledger(w->p);
    w->membranes0 = pl_vesicles_membranes(w->p);
    w->total0 = pl_vesicles_fatty_acid(w->p) + ledger.washed - ledger.renewed_monomer;
    /* from the first cycle on, a stopped seeded lineage keeps what its chemistry made and remakes nothing */
    if (w->o.stop && (status = pl_vesicles_set_inert(w->p, r, 1, d)) != PL_OK) return status;
    return step(w, d);
}
