/*
 * One protocell of Taneja and Higgs (Life 15, 2025) recorded for the viewer: an autocatalytic reaction inside a vesicle
 * of two leaflets swells it by osmosis, its membrane takes up lipid, it stretches into a capsule and divides when its
 * elastic energy drops; the viewer follows one daughter, the larger and the smaller in turn.
 *
 *   capsule [key=value ...]
 *     out=recordings/capsule-rflip-0.1 rflip=0.1 rmem=0.1 elip=1.01 (times C*) muw=0.015 equal=0 inactive=0
 *     duration=20000 every=20 dt=0.01 seed=5505
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "protolife/protolife.h"

typedef struct {
    const char *out;
    double rflip, rmem, elip, muw, duration, every, dt;
    int equal, inactive;
    unsigned long long seed;
} options;

static void parse(int argc, char **argv, options *o) {
    *o = (options){"recordings/capsule-rflip-0.1", 0.1, 0.1, 1.01, 0.015, 20000.0, 20.0, 0.01, 0, 0, 5505};
    for (int i = 1; i < argc; i++) {
        const char *eq = strchr(argv[i], '=');
        if (eq == NULL) continue;
        size_t n = (size_t)(eq - argv[i]);
        const char *v = eq + 1;
#define IS(key) (n == strlen(key) && strncmp(argv[i], key, n) == 0)
        if (IS("out")) o->out = v;
        else if (IS("rflip")) o->rflip = atof(v);
        else if (IS("rmem")) o->rmem = atof(v);
        else if (IS("elip")) o->elip = atof(v);
        else if (IS("muw")) o->muw = atof(v);
        else if (IS("duration")) o->duration = atof(v);
        else if (IS("every")) o->every = atof(v);
        else if (IS("dt")) o->dt = atof(v);
        else if (IS("equal")) o->equal = atoi(v);
        else if (IS("inactive")) o->inactive = atoi(v);
        else if (IS("seed")) o->seed = strtoull(v, NULL, 10);
#undef IS
    }
}

static const char *const VALUES[][3] = {
    {"volume", "volume", "V0"},
    {"lipid_area", "lipid area", "A0"},
    {"area", "area covered", "A0"},
    {"radius", "radius of the capsule", "R0"},
    {"length", "length of its cylinder", "R0"},
    {"reduced", "reduced volume", ""},
    {"energy", "elastic energy", "kappa"},
    {"food", "food C1 inside", "M"},
    {"catalyst", "catalyst C2 inside", "M"},
    {"waste", "waste W inside", "M"},
    {"lipid", "lipid in solution inside", "M"},
    {"divisions", "divisions so far", ""},
    {"bursts", "bursts so far", ""},
    {"event", "since the last frame: 1 a division, 2 a burst", ""},
    {"share", "volume share of the daughter followed, at the last division", ""}};
enum { VALUE_COUNT = sizeof VALUES / sizeof *VALUES };

static pl_status describe(const options *o, pl_scene **scene, pl_diagnostic *d) {
    char title[160], subtitle[320];
    snprintf(title, sizeof title, "%s protocell of Taneja and Higgs: flip-flop %g", o->inactive ? "An inactive" : "A", o->rflip);
    snprintf(subtitle, sizeof subtitle,
             "Lipid exchange %g, flip-flop %g, lipid outside %g C*, waste leaving at %g; %s; one daughter followed after "
             "each division, the larger and the smaller in turn",
             o->rmem, o->rflip, o->elip, o->muw, o->equal ? "division into equal halves at v = 1/2" : "division when the energy drops");
    pl_status status = pl_scene_create("capsule", title, subtitle,
                                       "protolife examples/capsule.c, after Taneja and Higgs, Life 15 (2025), Zenodo 15157714",
                                       scene, d);
    if (status == PL_OK) status = pl_scene_set_time_unit(*scene, "time units", d);
    for (int i = 0; i < VALUE_COUNT && status == PL_OK; i++) status = pl_scene_add_value(*scene, VALUES[i][0], VALUES[i][1], VALUES[i][2], d);
    const struct {
        const char *key, *label;
        double value;
    } constants[] = {{"rflip", "flip-flop rate", o->rflip}, {"rmem", "lipid exchange rate", o->rmem},
                     {"elip", "lipid outside, in C*", o->elip}, {"muw", "waste permeation", o->muw},
                     {"equal_division", "division into equal halves", o->equal}, {"inactive", "started without catalyst", o->inactive}};
    for (size_t i = 0; i < sizeof constants / sizeof *constants && status == PL_OK; i++)
        status = pl_scene_add_constant(*scene, constants[i].key, constants[i].label, "", constants[i].value, d);
    return status;
}

int main(int argc, char **argv) {
    options o;
    parse(argc, argv, &o);
    pl_network *network = NULL;
    pl_capsule *cell = NULL;
    pl_scene *scene = NULL;
    pl_diagnostic d = {0};
    pl_capsule_parameters k = pl_capsule_taneja_higgs();
    k.rflip = o.rflip, k.rmem = o.rmem, k.elip = o.elip * k.cstar, k.equal_division = o.equal;
    double mu[3] = {0.01, 0.0, o.muw}, outside[3] = {1.0, 0.0, 0.0};  /* food in, catalyst kept, waste out */
    pl_status status = pl_network_load(PROTOLIFE_DATA "/networks/taneja-higgs.net", &network, &d);
    if (status == PL_OK) status = pl_capsule_create(&k, network, mu, outside, o.seed, &cell, &d);
    if (status == PL_OK && o.inactive) pl_capsule_set_molar(cell, 0, 1.0, NULL), pl_capsule_set_molar(cell, 1, 0.0, NULL);
    if (status == PL_OK) status = pl_capsule_settle(cell, 1000.0, o.dt, &d);
    if (status == PL_OK) status = describe(&o, &scene, &d);
    if (status == PL_OK) status = pl_scene_open(scene, o.out, &d);
    long frames = lround(o.duration / o.every);
    int seen = 0;
    float row[VALUE_COUNT];
    for (long f = 0; status == PL_OK && f <= frames; f++) {
        if (f > 0) status = pl_capsule_advance(cell, o.every, o.dt, &d);
        pl_capsule_shape s = pl_capsule_get_shape(cell);
        const pl_capsule_event *e;
        int events = pl_capsule_events(cell, &e), divisions = 0, bursts = 0, event = 0;
        double share = NAN;
        for (int i = 0; i < events; i++) {
            if (e[i].kind == PL_CAPSULE_DIVISION) divisions++, share = e[i].share;
            else bursts++;
            if (i >= seen) event = event == 1 ? 1 : (e[i].kind == PL_CAPSULE_DIVISION ? 1 : 2);
        }
        seen = events;
        const double v[VALUE_COUNT] = {s.volume, s.lipid_area, s.area, s.radius, s.length, s.reduced, s.energy,
                                       pl_capsule_molar(cell, 0), pl_capsule_molar(cell, 1), pl_capsule_molar(cell, 2),
                                       pl_capsule_lipid(cell), divisions, bursts, event, share};
        for (int i = 0; i < VALUE_COUNT; i++) row[i] = (float)v[i];
        pl_status written = pl_scene_frame(scene, pl_capsule_time(cell), row, NULL, &d);
        if (status == PL_OK) status = written;
    }
    if (status != PL_OK) fprintf(stderr, "%s at t = %g: %s\n", pl_status_name(status), d.time, d.message);
    else printf("%s: %ld frames\n", o.out, pl_scene_frames(scene));
    pl_scene_free(scene);
    pl_capsule_free(cell);
    pl_network_free(network);
    return status == PL_OK || status == PL_ERROR_MODEL ? 0 : 1;  /* a cell that ends as the authors' do is a result */
}
