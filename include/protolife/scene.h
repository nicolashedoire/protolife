#ifndef PROTOLIFE_SCENE_H
#define PROTOLIFE_SCENE_H

#include "protolife/diagnostic.h"

/*
 * A recording of a simulation for a viewer (the Godot viewer in viewer/): a directory holding scene.json, which says
 * what the world is and what each number means, and frames.bin, frames of a fixed size appended as the simulation
 * runs, so that a viewer can follow a run live. A frame is the time (float64, s), then the values (float32 each, in the
 * manifest's order), then the fields (float32, field after field, each over the grid's cells, row-major). The format,
 * "protolife-scene-1", is Genesis's "genesis-scene-1" under another name.
 *
 * Build the manifest (kind, titles, values, fields, constants), open it in a directory, then append frames.
 */

typedef struct pl_scene pl_scene;

/* A manifest of this kind ("lineage", "vesicles", "vessel", "world": what the viewer draws), with a title, a subtitle
 * and the sources it rests on. */
pl_status pl_scene_create(const char *kind, const char *title, const char *subtitle, const char *source,
                          pl_scene **scene, pl_diagnostic *diagnostic);
/* A square grid of n x n cells over `box` metres, for fields (1 and 0 for a recording without fields: the default). */
pl_status pl_scene_set_grid(pl_scene *scene, int n, double box, pl_diagnostic *diagnostic);
/* A number of each frame: a stable key ("r0_radius"), a label for the screen and a unit. */
pl_status pl_scene_add_value(pl_scene *scene, const char *key, const char *label, const char *unit,
                             pl_diagnostic *diagnostic);
/* A field of each frame, over the grid. */
pl_status pl_scene_add_field(pl_scene *scene, const char *key, const char *label, const char *unit,
                             pl_diagnostic *diagnostic);
/* A number of the whole recording. */
pl_status pl_scene_add_constant(pl_scene *scene, const char *key, const char *label, const char *unit, double value,
                                pl_diagnostic *diagnostic);
/* Makes the directory and its parents, writes scene.json and opens frames.bin afresh. The manifest is then fixed. */
pl_status pl_scene_open(pl_scene *scene, const char *directory, pl_diagnostic *diagnostic);
/* Appends one frame (values: one per value; fields: cells per field, NULL without fields) and flushes it. */
pl_status pl_scene_frame(pl_scene *scene, double time, const float *values, const float *fields,
                         pl_diagnostic *diagnostic);
int pl_scene_value_count(const pl_scene *scene);
long pl_scene_frames(const pl_scene *scene);  /* frames written */
/* Closes the recording and frees the manifest. */
void pl_scene_free(pl_scene *scene);

#endif
