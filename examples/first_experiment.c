/*
 * A first experiment with protolife: does a chemistry inside make a vesicle grow and divide? Two protocells of Taneja
 * and Higgs (Life 15, 2025) under the same conditions: one holds an autocatalytic reaction, its control holds the
 * same food but no catalyst. The program prints them side by side and records both for the viewer. It runs in a few
 * seconds; change the settings below and run it again.
 *
 *   first_experiment                      (writes recordings/first-active and recordings/first-control)
 *
 * Read docs/getting-started.md for a tour of this file.
 */
#include <stdio.h>

#include "protolife/protolife.h"

/* ---- the settings: change them ---- */
static const double LIPID_OUTSIDE = 1.01;    /* lipid in the solution outside, in units of the balance C*: above 1 the
                                                membrane can grow */
static const double FLIP_FLOP = 0.1;         /* lipid flipping between the two leaflets: above 0.18 no division */
static const double WASTE_LEAVES = 0.015;    /* how fast the waste leaves: at or below the food's 0.01, no steady size */
static const double DURATION = 20000.0;      /* time units */
static const double REPORT_EVERY = 2000.0;   /* a line of the table */

/* One protocell: the standard set of Taneja and Higgs with the settings above. A control starts without catalyst. */
static pl_capsule *protocell(const pl_network *chemistry, int control, pl_diagnostic *d) {
    pl_capsule_parameters k = pl_capsule_taneja_higgs();
    k.elip = LIPID_OUTSIDE * k.cstar, k.rflip = FLIP_FLOP;
    /* per species of the network (food C1, catalyst C2, waste W): how fast it crosses, and what is outside */
    double crossing[3] = {0.01, 0.0, WASTE_LEAVES}, outside[3] = {1.0, 0.0, 0.0};
    pl_capsule *cell = NULL;
    if (pl_capsule_create(&k, chemistry, crossing, outside, 5505, &cell, d) != PL_OK) return NULL;
    if (control) {  /* the same cell with food only: no catalyst, so no reaction */
        pl_capsule_set_molar(cell, pl_network_find(chemistry, "C1"), 1.0, d);
        pl_capsule_set_molar(cell, pl_network_find(chemistry, "C2"), 0.0, d);
    }
    pl_capsule_settle(cell, 1000.0, 0.01, d);  /* the chemistry settles before the membrane moves */
    return cell;
}

static int divisions(const pl_capsule *cell) {
    const pl_capsule_event *events;
    int n = pl_capsule_events(cell, &events), count = 0;
    for (int i = 0; i < n; i++) count += events[i].kind == PL_CAPSULE_DIVISION;
    return count;
}

/* A recording for the viewer: what each number of a frame means. */
static pl_scene *recording(const char *directory, const char *title, pl_diagnostic *d) {
    static const char *const keys[][3] = {{"volume", "volume", "V0"}, {"lipid_area", "lipid area", "A0"},
                                          {"area", "area covered", "A0"}, {"radius", "radius", "R0"},
                                          {"length", "length of its cylinder", "R0"}, {"reduced", "reduced volume", ""},
                                          {"energy", "elastic energy", "kappa"}, {"food", "food inside", "M"},
                                          {"catalyst", "catalyst inside", "M"}, {"waste", "waste inside", "M"},
                                          {"lipid", "lipid inside", "M"}, {"divisions", "divisions so far", ""},
                                          {"bursts", "bursts so far", ""}, {"event", "1 a division, 2 a burst", ""},
                                          {"share", "share of the daughter followed", ""}};
    pl_scene *scene = NULL;
    if (pl_scene_create("capsule", title, "protolife examples/first_experiment.c", "Taneja and Higgs, Life 15 (2025)",
                        &scene, d) != PL_OK)
        return NULL;
    pl_scene_set_time_unit(scene, "time units", d);
    for (size_t i = 0; i < sizeof keys / sizeof *keys; i++) pl_scene_add_value(scene, keys[i][0], keys[i][1], keys[i][2], d);
    if (pl_scene_open(scene, directory, d) != PL_OK) return pl_scene_free(scene), NULL;
    return scene;
}

static void frame(pl_scene *scene, const pl_capsule *cell, int *seen) {
    pl_capsule_shape s = pl_capsule_get_shape(cell);
    const pl_capsule_event *events;
    int n = pl_capsule_events(cell, &events), bursts = 0, event = 0;
    double share = 0.0;
    for (int i = 0; i < n; i++) {
        bursts += events[i].kind == PL_CAPSULE_BURST;
        if (events[i].kind == PL_CAPSULE_DIVISION) share = events[i].share;
        if (i >= *seen) event = event == 1 || events[i].kind == PL_CAPSULE_DIVISION ? 1 : 2;
    }
    *seen = n;
    float row[15] = {(float)s.volume, (float)s.lipid_area, (float)s.area, (float)s.radius, (float)s.length,
                     (float)s.reduced, (float)s.energy, (float)pl_capsule_molar(cell, 0), (float)pl_capsule_molar(cell, 1),
                     (float)pl_capsule_molar(cell, 2), (float)pl_capsule_lipid(cell), (float)divisions(cell),
                     (float)bursts, (float)event, (float)share};
    pl_scene_frame(scene, pl_capsule_time(cell), row, NULL, NULL);
}

int main(void) {
    pl_diagnostic d = {0};
    pl_network *chemistry = NULL;
    if (pl_network_load(PROTOLIFE_DATA "/networks/taneja-higgs.net", &chemistry, &d) != PL_OK) {
        fprintf(stderr, "%s\n", d.message);
        return 1;
    }
    pl_capsule *active = protocell(chemistry, 0, &d), *control = protocell(chemistry, 1, &d);
    pl_scene *scenes[2] = {recording("recordings/first-active", "First experiment: the active protocell", &d),
                           recording("recordings/first-control", "First experiment: the control without catalyst", &d)};
    if (active == NULL || control == NULL || scenes[0] == NULL || scenes[1] == NULL) {
        fprintf(stderr, "%s\n", d.message);
        return 1;
    }
    printf("%10s | %-34s | %-34s\n", "time", "active: volume, lipid area, divisions", "control: volume, lipid area, divisions");
    pl_capsule *cells[2] = {active, control};
    int seen[2] = {0, 0}, ended[2] = {0, 0};
    for (double t = 0.0; t <= DURATION + 1e-9; t += 20.0) {
        for (int c = 0; c < 2; c++) {
            if (t > 0.0 && !ended[c] && pl_capsule_advance(cells[c], 20.0, 0.01, &d) != PL_OK) {
                ended[c] = 1;  /* the model ends there, as the authors' program does: a thin tube, or unbounded growth */
                printf("%10.0f | %s cell: %s\n", t, c == 0 ? "the active" : "the control", d.message);
            }
            if (!ended[c]) frame(scenes[c], cells[c], &seen[c]);
        }
        if (t - REPORT_EVERY * (long)(t / REPORT_EVERY) < 1e-9) {
            pl_capsule_shape a = pl_capsule_get_shape(active), b = pl_capsule_get_shape(control);
            printf("%10.0f | %10.3f %10.3f %10d     | %10.3f %10.3f %10d\n", t, a.volume, a.lipid_area, divisions(active),
                   b.volume, b.lipid_area, divisions(control));
        }
    }
    printf("\nRecorded: recordings/first-active and recordings/first-control. Watch them: godot --path viewer\n");
    pl_scene_free(scenes[0]), pl_scene_free(scenes[1]);
    pl_capsule_free(active), pl_capsule_free(control);
    pl_network_free(chemistry);
    return 0;
}
