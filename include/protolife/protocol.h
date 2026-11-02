#ifndef PROTOLIFE_PROTOCOL_H
#define PROTOLIFE_PROTOCOL_H

#include "protolife/vesicles.h"

/*
 * Laboratory protocols applied to a vesicle population as time passes, after Genesis world/vesicle_retention.c.
 *
 * Cycles of growth and division (after Zhu and Szostak, J. Am. Chem. Soc. 131, 5705, 2009, who dilute 1:10 before each
 * cycle): at the start of each cycle a share of the vessel is replaced by fresh medium; then micelles are fed in doses
 * that together multiply the membranes by 1 + `micelles`, the first 10 s into the cycle and the others `interval` s
 * apart; and the vessel is sheared from two thirds of the cycle for 30 s, or from the first cycle on when `stirred`.
 */
typedef struct pl_cycles {
    int count;        /* cycles; 0: none */
    double start;     /* s: when the first one begins */
    double period;    /* s */
    double dilution;  /* share of the vessel replaced at the start of each cycle */
    double micelles;  /* equivalents of the membranes fed per cycle, all doses together */
    int feeds;        /* doses */
    double interval;  /* s between doses */
    double shear;     /* 1/s */
    int stirred;      /* 1: sheared from the start of the first cycle on */
} pl_cycles;

/* Applies what the cycles do while time goes from `before` to `after` (s): call it after each step. */
void pl_cycles_apply(pl_vesicles *population, const pl_cycles *cycles, double before, double after);
/* The share of the vessel kept through the dilutions made by `time`. */
double pl_cycles_kept(const pl_cycles *cycles, double time);

#endif
