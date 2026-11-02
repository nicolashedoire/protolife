/*
 * Living vesicles through cycles of feeding and shear, recorded for the viewer: step 6 of Genesis's goal (entries 83
 * to 85). A network of Origin's sugar chemistry (an E0 core) runs inside 500 nm oleate vesicles and in their medium, in
 * a dialysis bag renewed with formaldehyde; most of the fatty acid starts in vesicles holding the core's state, the
 * rest in vesicles of fresh medium. After an hour of installation, unrecorded, each cycle halves the bag, feeds one
 * equivalent of micelles and shears the vessel. With living=0 no vesicle reacts: the control.
 *
 *   reproduction [key=value ...]
 *     network=data/networks/e0-random2.net out=recordings/step-6-random2 name=random2 living=1 seed=core|fresh stop=0
 *     vesicles=0 (a finite bag of that many seeded vesicles, its dilutions drawn at random) draw=1
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
    *o = defaults();
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
        else if (IS("seed")) o->twins = strcmp(v, "fresh") == 0;
        else if (IS("stop")) o->stop = (int)x;
        else if (IS("vesicles")) o->vesicles = x;
        else if (IS("draw")) o->draw = strtoull(v, NULL, 10);
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
    dismantle(&w);
    return status == PL_OK ? 0 : 1;
}
