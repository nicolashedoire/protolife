/*
 * Living vesicles through cycles of feeding and shear, recorded for the viewer: step 6 of Genesis's goal (entries 83
 * to 85). A network of Origin's sugar chemistry (an E0 core) runs inside 500 nm oleate vesicles and in their medium, in
 * a dialysis bag renewed with formaldehyde; most of the fatty acid starts in vesicles holding the core's state, the
 * rest in vesicles of fresh medium. After an hour of installation, unrecorded, each cycle halves the bag, feeds one
 * equivalent of micelles and shears the vessel. With living=0 no vesicle reacts: the control.
 *
 *   reproduction [key=value ...]
 *     network=data/networks/e0-random2.net out=recordings/step-6-random2 name=random2 living=1
 *     radius=500 (nm) acid=3 (mM) seeded=0.9 small=2.9e-8 (cm/s) food=0.5 (M) renewal=0.05 (1/s)
 *     cycles=4 period=900 micelles=1 feeds=1 interval=60 shear=20 stirred=1 start=3600 dilution=0.5
 *     frames=361 every=10 dt=0.5
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "reproduction.h"

static void parse(int argc, char **argv, options *o) {
    *o = (options){"data/networks/e0-random2.net", "recordings/step-6-random2", "random2", 1, 500e-9, 3e-3, 0.9,
                   2.9e-10, 0.5, 5e-2, 361, 10.0, 0.5, {4, 3600.0, 900.0, 0.5, 1.0, 1, 60.0, 20.0, 1}};
    for (int i = 1; i < argc; i++) {
        const char *eq = strchr(argv[i], '=');
        if (eq == NULL) continue;
        size_t n = (size_t)(eq - argv[i]);
        const char *v = eq + 1;
        double x = atof(v);
#define IS(key) (n == strlen(key) && strncmp(argv[i], key, n) == 0)
        if (IS("network")) o->network = v;
        else if (IS("out")) o->out = v;
        else if (IS("name")) o->name = v;
        else if (IS("living")) o->living = (int)x;
        else if (IS("radius")) o->radius = x * 1e-9;
        else if (IS("acid")) o->acid = x * 1e-3;
        else if (IS("seeded")) o->seeded = x;
        else if (IS("small")) o->small = x * 1e-2;
        else if (IS("food")) o->food = x;
        else if (IS("renewal")) o->renewal = x;
        else if (IS("frames")) o->frames = (long)x;
        else if (IS("every")) o->every = x;
        else if (IS("dt")) o->dt = x;
        else if (IS("cycles")) o->cycles.count = (int)x;
        else if (IS("period")) o->cycles.period = x;
        else if (IS("micelles")) o->cycles.micelles = x;
        else if (IS("feeds")) o->cycles.feeds = (int)x;
        else if (IS("interval")) o->cycles.interval = x;
        else if (IS("shear")) o->cycles.shear = x;
        else if (IS("stirred")) o->cycles.stirred = (int)x;
        else if (IS("start")) o->cycles.start = x;
        else if (IS("dilution")) o->cycles.dilution = x;
#undef IS
    }
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
    filled(w, core, volume, inside);
    pl_status status = pl_vesicles_add(w->p, o->acid * o->seeded * AVOGADRO / molecules, molecules, inside, SEEDED, NULL, d);
    filled(w, fresh, volume, inside);
    if (status == PL_OK)
        status = pl_vesicles_add(w->p, o->acid * (1.0 - o->seeded) * AVOGADRO / molecules, molecules, inside, FRESH, NULL, d);
    free(core), free(fresh), free(inside);
    return status;
}

static pl_status build(world *w, pl_diagnostic *d) {
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

pl_status step(world *w, pl_diagnostic *d) {
    double before = pl_vesicles_time(w->p);
    pl_status status = pl_vesicles_advance(w->p, w->o.dt, d);
    if (status == PL_OK) pl_cycles_apply(w->p, &w->o.cycles, before, pl_vesicles_time(w->p));
    return status;
}

/* The installation, unrecorded, up to a step before the first cycle; there, the reference; then that step. */
static pl_status install(world *w, pl_diagnostic *d) {
    pl_status status = PL_OK;
    while (status == PL_OK && pl_vesicles_time(w->p) + w->o.dt < w->o.cycles.start - w->o.dt / 2.0) status = step(w, d);
    int r = pl_vesicles_find(w->p, SEEDED);
    if (status == PL_OK && r < 0) return PL_ERROR_MODEL;  /* the seeded vesicles came apart */
    if (status != PL_OK) return status;
    w->count0 = pl_vesicles_record(w->p, r).count, w->held0 = sugar_carbon(w, r);
    pl_vesicle_ledger ledger = pl_vesicles_ledger(w->p);
    w->membranes0 = pl_vesicles_membranes(w->p);
    w->total0 = pl_vesicles_fatty_acid(w->p) + ledger.washed - ledger.renewed_monomer;
    return step(w, d);
}

int main(int argc, char **argv) {
    world w = {0};
    parse(argc, argv, &w.o);
    pl_diagnostic d = {0};
    pl_status status = build(&w, &d);
    if (status == PL_OK) status = install(&w, &d);
    pl_scene *scene = NULL;
    if (status == PL_OK) status = describe(&w, &scene, &d);
    if (status == PL_OK) status = pl_scene_open(scene, w.o.out, &d);
    if (status == PL_OK) status = record(&w, scene, &d);
    if (status != PL_OK) fprintf(stderr, "%s at t = %.6g s: %s\n", pl_status_name(status), d.time, d.message);
    else printf("%s: %ld frames: %s\n", w.o.out, pl_scene_frames(scene), w.o.living ? "living vesicles" : "the control");
    const pl_balance *balance = NULL;
    int balances = w.p != NULL ? pl_vesicles_balances(w.p, &balance) : 0;
    for (int b = 0; b < balances; b++)  /* the balances checked at every step, and how far they drifted in all */
        printf("  balance of %s: %.6g mol/L, drift %.3g (%.3g of it)\n", balance[b].name, balance[b].value,
               balance[b].drift, balance[b].value != 0.0 ? fabs(balance[b].drift / balance[b].value) : 0.0);
    pl_scene_free(scene);
    pl_vesicles_free(w.p);
    pl_network_free(w.network);
    return status == PL_OK ? 0 : 1;
}
