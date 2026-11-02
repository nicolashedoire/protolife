/* The recording of the reproduction example for the viewer, kind "lineage" (reproduction.h), after Genesis
 * world/scene_reproduction.c. Three groups, by what they are rather than by record: r0_ the seeded vesicles and their
 * descendants, r1_ the fresh vesicles of the start, r2_ every vesicle born of micelles. Counts are per seeded vesicle
 * of the first cycle, corrected for the dilutions. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "reproduction.h"

static const char *const GROUP[GROUPS] = {"the lineage", "the fresh vesicles of the start", "the vesicles born of micelles"};
static const char *const WHAT[PER_GROUP][3] = {
    {"present", "exists", ""},
    {"count", "vesicles per seeded vesicle", ""},
    {"generation", "generation", ""},
    {"origin", "origin (0: the start; 1: the micelles)", ""},
    {"radius", "radius of the sphere of its area", "nm"},
    {"reduced", "reduced volume", ""},
    {"tension", "tension", "mN/m"},
    {"cargo_share", "sugar carbon per vesicle, against a seeded vesicle at the first cycle", ""},
    {"cargo", "sugar carbon inside", "M"},
    {"share", "share of the membranes", ""}};
static const char *const EXTRA[][3] = {
    {"shear", "shear", "1/s"},
    {"coat", "micelles on the membranes", "mM"},
    {"intermediate", "intermediate going into the membranes", "mM"},
    {"forming", "intermediate making new vesicles", "mM"},
    {"born", "vesicles born of micelles, per seeded vesicle", ""},
    {"divided", "vesicles of the lineage made by division, per seeded vesicle", ""},
    {"fed", "fatty acid added, in equivalents of the membranes at the first cycle", ""},
    {"monomer", "free fatty acid", "mM"},
    {"medium_sugar", "sugar carbon in the medium", "M"}};

static pl_status values_and_constants(const world *w, pl_scene *s, pl_diagnostic *d) {
    char key[48], label[160];
    pl_status status = PL_OK;
    for (int g = 0; g < GROUPS && status == PL_OK; g++)
        for (int q = 0; q < PER_GROUP && status == PL_OK; q++) {
            snprintf(key, sizeof key, "r%d_%s", g, WHAT[q][0]);
            snprintf(label, sizeof label, "%s: %s", GROUP[g], WHAT[q][1]);
            status = pl_scene_add_value(s, key, label, WHAT[q][2], d);
        }
    for (size_t i = 0; i < sizeof EXTRA / sizeof *EXTRA && status == PL_OK; i++)
        status = pl_scene_add_value(s, EXTRA[i][0], EXTRA[i][1], EXTRA[i][2], d);
    const pl_cycles *c = &w->o.cycles;
    pl_vesicle_constants k = pl_vesicle_oleate();
    const struct {
        const char *key, *label, *unit;
        double value;
    } constants[] = {
        {"records", "groups of vesicles", "", GROUPS},
        {"inside_sugars", "the inside is made by a chemistry (Origin's sugars)", "", 1.0},
        {"living", "the chemistry runs in the vesicles", "", w->o.living},
        {"micelles", "fed per cycle (equivalents of the membranes)", "", c->micelles},
        {"feeds", "doses", "", c->feeds},
        {"interval", "between two doses", "s", c->interval},
        {"first_feed", "first dose", "s", 10.0},
        {"cycles", "cycles", "", c->count},
        {"period", "length of a cycle", "s", c->period},
        {"stirred", "sheared all along", "", c->stirred},
        {"shear_at", "shear, in each cycle, at", "s", c->stirred ? 0.0 : 2.0 * c->period / 3.0},
        {"shear_for", "shear for", "s", c->stirred ? c->period : 30.0},
        {"shear_rate", "shear applied", "1/s", c->shear},
        {"dilution", "share of the bag renewed at the start of each cycle", "", c->dilution},
        {"installation", "installation before the first cycle, not recorded", "s", c->start},
        {"critical_shear", "critical shear", "1/s", k.critical_shear},
        {"division_reduced", "reduced volume below which a sheared vesicle divides", "", k.division_reduced},
        {"least_radius", "smallest daughter", "nm", k.least_radius * 1e9},
        {"born_radius", "radius of the vesicles born of micelles", "nm", k.born_radius * 1e9},
        {"rupture", "rupture tension", "mN/m", k.rupture * 1000.0}};
    for (size_t i = 0; i < sizeof constants / sizeof *constants && status == PL_OK; i++)
        status = pl_scene_add_constant(s, constants[i].key, constants[i].label, constants[i].unit, constants[i].value, d);
    return status;
}

pl_status describe(const world *w, pl_scene **scene, pl_diagnostic *d) {
    char title[160], subtitle[400];
    const pl_cycles *c = &w->o.cycles;
    snprintf(title, sizeof title, w->o.living ? "Living vesicles: core %s" : "The control without chemistry: core %s",
             w->o.name);
    snprintf(subtitle, sizeof subtitle,
             "Oleate vesicles of %.3g nm in a dialysis bag, formaldehyde %.3g M; %d cycles of %.0f min (after %.0f min "
             "of installation); %s; shear %s",
             w->o.radius * 1e9, w->o.food, c->count, c->period / 60.0, c->start / 60.0,
             c->feeds > 1 ? "micelles in slow doses" : "one dose of micelles at once",
             c->stirred ? "all along" : "at two thirds of each cycle");
    pl_status status = pl_scene_create("lineage", title, subtitle,
                                       "protolife examples/reproduction.c, after Genesis entries 83 to 85; Zhu and "
                                       "Szostak 2009; Chen and Szostak 2004",
                                       scene, d);
    return status == PL_OK ? values_and_constants(w, *scene, d) : status;
}

static int group_of(int label) { return label == SEEDED ? LINEAGE : label == FRESH ? FRESH_GROUP : label == 0 ? BORN : -1; }

typedef struct {
    double count[GROUPS], area[GROUPS], volume[GROUPS], sphere[GROUPS], tension[GROUPS], held[GROUPS], acid[GROUPS];
    int generation[GROUPS], present[GROUPS];
} sums;

static void add_up(const world *w, sums *t) {
    memset(t, 0, sizeof *t);
    for (int r = 0; r < pl_vesicles_records(w->p); r++) {
        pl_vesicle v = pl_vesicles_record(w->p, r);
        int g = group_of(v.label);
        if (g < 0) continue;
        pl_vesicle_shape s = pl_vesicles_shape(w->p, r);
        t->present[g] = 1;
        t->count[g] += v.count, t->area[g] += v.count * s.area, t->volume[g] += v.count * s.volume;
        t->sphere[g] += v.count * 4.0 / 3.0 * PI * pow(s.radius, 3.0), t->tension[g] += v.count * s.tension;
        t->held[g] += v.count * sugar_carbon(w, r);
        t->acid[g] += v.count * v.molecules / AVOGADRO;
        if (v.generation > t->generation[g]) t->generation[g] = v.generation;
    }
}

static void values(const world *w, float *out) {
    double per = w->count0 * pl_cycles_kept(&w->o.cycles, pl_vesicles_time(w->p)), membranes = pl_vesicles_membranes(w->p);
    double coat = 0.0, intermediate = 0.0;
    for (int r = 0; r < pl_vesicles_records(w->p); r++) {
        pl_vesicle v = pl_vesicles_record(w->p, r);
        coat += v.coat, intermediate += v.intermediate;
    }
    sums t;
    add_up(w, &t);
    int at = 0;
    for (int g = 0; g < GROUPS; g++) {
        out[at++] = (float)t.present[g];
        if (!t.present[g]) {
            for (int q = 1; q < PER_GROUP; q++) out[at++] = NAN;
            continue;
        }
        out[at++] = (float)(t.count[g] / per);
        out[at++] = (float)t.generation[g];
        out[at++] = g == BORN ? 1.0f : 0.0f;
        out[at++] = (float)(sqrt(t.area[g] / t.count[g] / (4.0 * PI)) * 1e9);
        out[at++] = (float)(t.volume[g] / t.sphere[g]);
        out[at++] = (float)(t.tension[g] / t.count[g] * 1000.0);
        out[at++] = (float)(w->held0 > 0.0 ? t.held[g] / t.count[g] / w->held0 : NAN);
        out[at++] = (float)(t.held[g] / (t.volume[g] * 1000.0));
        out[at++] = (float)(t.acid[g] / membranes);
    }
    pl_vesicle_ledger ledger = pl_vesicles_ledger(w->p);
    out[at++] = (float)pl_vesicles_shear(w->p);
    out[at++] = (float)(coat * 1000.0);
    out[at++] = (float)(intermediate * 1000.0);
    out[at++] = (float)(pl_vesicles_forming(w->p) * 1000.0);
    out[at++] = (float)(t.count[BORN] / per);
    out[at++] = (float)(t.present[LINEAGE] ? t.count[LINEAGE] / per - 1.0 : NAN);
    out[at++] = (float)((pl_vesicles_fatty_acid(w->p) + ledger.washed - ledger.renewed_monomer - w->total0) / w->membranes0);
    out[at++] = (float)(1000.0 * pl_vesicles_monomer(w->p) / pl_vesicles_medium(w->p, NULL));
    out[at++] = (float)sugar_carbon(w, -1);
}

pl_status record(world *w, pl_scene *scene, pl_diagnostic *d) {
    float row[GROUPS * PER_GROUP + sizeof EXTRA / sizeof *EXTRA];
    long per_frame = (long)llround(w->o.every / w->o.dt);
    pl_status status = PL_OK;
    for (long frame = 0; status == PL_OK && frame < w->o.frames; frame++) {
        for (long s = 0; status == PL_OK && frame > 0 && s < (per_frame > 0 ? per_frame : 1); s++) status = step(w, d);
        if (status != PL_OK) break;
        values(w, row);
        status = pl_scene_frame(scene, pl_vesicles_time(w->p) - w->o.cycles.start, row, NULL, d);
    }
    return status;
}
