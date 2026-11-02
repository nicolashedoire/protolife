/* Laboratory protocols on vesicle populations (protolife/protocol.h). */
#include "protolife/protocol.h"

#include <math.h>
#include <stddef.h>

static int passes(double moment, double before, double after) { return before < moment && after >= moment; }

void pl_cycles_apply(pl_vesicles *p, const pl_cycles *c, double before, double after) {
    if (p == NULL || c == NULL) return;
    for (int k = 0; k < c->count; k++) {
        double start = c->start + k * c->period, shear = c->stirred ? c->start : start + 2.0 * c->period / 3.0;
        if (c->dilution > 0.0 && passes(start, before, after)) pl_vesicles_dilute(p, c->dilution);
        for (int i = 0; i < c->feeds; i++)  /* doses that together multiply the membranes by 1 + micelles */
            if (passes(start + 10.0 + i * c->interval, before, after))
                pl_vesicles_feed_micelles(p, (pow(1.0 + c->micelles, 1.0 / c->feeds) - 1.0) * pl_vesicles_membranes(p));
        if (passes(shear, before, after)) pl_vesicles_set_shear(p, c->shear, NULL);
        if (!c->stirred && passes(shear + 30.0, before, after)) pl_vesicles_set_shear(p, 0.0, NULL);
    }
}

double pl_cycles_kept(const pl_cycles *c, double time) {
    int done = 0;
    for (int k = 0; c != NULL && k < c->count; k++) done += time >= c->start + k * c->period;
    return c != NULL ? pow(1.0 - c->dilution, done) : 1.0;
}
