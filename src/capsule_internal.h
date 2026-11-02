#ifndef PROTOLIFE_CAPSULE_INTERNAL_H
#define PROTOLIFE_CAPSULE_INTERNAL_H

#include "protolife/capsule.h"

/* The form a volume takes in a lipid area (both relative to the sphere of radius R0), after Taneja and Higgs's
 * growADE.c: a sphere when the volume fills the area, else a capsule found by bisection on its radius. */
typedef struct {
    double radius, length;  /* R/R0, L/R0 */
    double actual;          /* the area the membrane covers: the lipid area, or the swollen sphere's */
    int swollen;
} capsule_form;

capsule_form pl_capsule_form(double volume, double lipid_area);

/* The elastic energy of a cell of this volume and leaflet areas (relative), in units of kappa's energy. DA0 is the
 * area difference of the sphere of radius R0 over its area, 4 d / R0. */
double pl_capsule_energy(const pl_capsule_parameters *k, double DA0, double volume, double outer, double inner);

#endif
