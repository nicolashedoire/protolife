/* Step 8 of Genesis's goal, shortened to one cycle (Genesis entries 90 and 91, through examples/reproduction_world.c):
 * half the fatty acid in vesicles holding the E0 carbons core, half in vesicles of fresh medium; half an hour of
 * installation, then a cycle of dilution, micelles and shear. The argument names the arm:
 *   living   the core runs: its lineage's share of the membranes rises (by 0.118 when written);
 *   twins    the seeded vesicles hold fresh medium too: the share stays one half exactly, all along;
 *   stopped  the core stops at the cycle: the share stays put (it moved by -0.0007 when written).
 * One arm per run, so that ctest runs them side by side. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "check.h"
#include "reproduction.h"

static double share(const world *w) {
    double seeded = fatty_acid_of(w, SEEDED), fresh = fatty_acid_of(w, FRESH);
    return seeded / (seeded + fresh);
}

/* The share at the cycle's start and at its end, and its largest gap from one half during the cycle. */
static void run(int twins, int stop, double *start, double *end, double *gap) {
    world w = {.o = defaults()};
    w.o.network = PROTOLIFE_DATA "/networks/e0-carbons.net";
    w.o.seeded = 0.5, w.o.twins = twins, w.o.stop = stop;
    w.o.cycles.start = 1800.0, w.o.cycles.count = 1, w.o.cycles.period = 900.0;
    pl_diagnostic d = {0};
    pl_status status = build(&w, &d);
    if (status == PL_OK) status = install(&w, &d);
    *start = share(&w), *gap = fabs(*start - 0.5);
    double finish = w.o.cycles.start + w.o.cycles.period - w.o.dt / 2.0;
    while (status == PL_OK && pl_vesicles_time(w.p) < finish) {
        status = step(&w, &d);
        *gap = fmax(*gap, fabs(share(&w) - 0.5));
    }
    if (status != PL_OK) fprintf(stderr, "%s at t = %.6g s: %s\n", pl_status_name(status), d.time, d.message);
    CHECK(status == PL_OK);
    *end = share(&w);
    printf("share %.6f at the cycle, %.6f after it; largest gap from one half %.3g\n", *start, *end, *gap);
    dismantle(&w);
}

int main(int argc, char **argv) {
    const char *arm = argc > 1 ? argv[1] : "";
    double start, end, gap;
    if (strcmp(arm, "living") == 0) {
        run(0, 0, &start, &end, &gap);
        CHECK(start > 0.6 && end - start > 0.05);  /* the living lineage took more than the fresh one, and takes more */
    } else if (strcmp(arm, "twins") == 0) {
        run(1, 0, &start, &end, &gap);
        CHECK(gap < 1e-12);  /* the simulation favours no record: twins stay twins */
    } else if (strcmp(arm, "stopped") == 0) {
        run(0, 1, &start, &end, &gap);
        CHECK(start > 0.6 && fabs(end - start) < 0.01);  /* what it made before stays; it gains nothing more */
    } else {
        fprintf(stderr, "test_selection living|twins|stopped\n");
        return 2;
    }
    return check_report(arm);
}
