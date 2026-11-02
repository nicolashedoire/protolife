#ifndef PROTOLIFE_EXAMPLE_REPRODUCTION_H
#define PROTOLIFE_EXAMPLE_REPRODUCTION_H

/* The world of the reproduction example (reproduction.c) and its recording (reproduction_scene.c). */
#include "protolife/protolife.h"

#define PI 3.14159265358979323846
#define AVOGADRO 6.02214076e23
#define BUFFER_MOLAR 0.2      /* bicine */
#define BUFFER_OSMOTIC 1.635  /* osmoles per mole of bicine at pH 8.5: the zwitterion and its sodium salt */
#define EXCESS 0.27           /* the relaxed vesicles' area above the sphere of their volume, a share of it */
enum { SEEDED = 1, FRESH = 2 };  /* labels; newborns have 0 */
enum { LINEAGE, FRESH_GROUP, BORN, GROUPS };
#define PER_GROUP 10

typedef struct {
    const char *network, *out, *name;
    int living;
    double radius, acid, seeded, small, food, renewal;
    long frames;
    double every, dt;
    pl_cycles cycles;
    int twins;  /* 1: the seeded vesicles hold fresh medium too, the twins of the fresh ones (a neutral control) */
    int stop;   /* 1: the seeded lineage's chemistry stops at the first cycle: it keeps what it made, remakes nothing */
    double vesicles;  /* a finite bag holding this many seeded vesicles at the start, its dilutions drawn; 0: mean field */
    unsigned long long draw;  /* the seed of those draws */
} options;

typedef struct {
    options o;
    pl_network *network;
    pl_vesicles *p;
    int buffer, food;
    double count0, held0, membranes0, total0;  /* a seeded vesicle at the first cycle, and the fatty acid then */
} world;

/* The manifest of the recording; then the frames, running the world and its cycles. */
pl_status describe(const world *w, pl_scene **scene, pl_diagnostic *diagnostic);
pl_status record(world *w, pl_scene *scene, pl_diagnostic *diagnostic);
/* One step and the cycles' events; carbon in sugars held by one vesicle of a record (mol), or in the medium (mol/L)
 * for record -1. */
pl_status step(world *w, pl_diagnostic *diagnostic);
double sugar_carbon(const world *w, int record);

#endif
