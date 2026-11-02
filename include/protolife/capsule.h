#ifndef PROTOLIFE_CAPSULE_H
#define PROTOLIFE_CAPSULE_H

#include <stdint.h>

#include "protolife/diagnostic.h"
#include "protolife/network.h"

/*
 * One protocell after Taneja and Higgs, "Protocell dynamics: modelling growth and division of lipid vesicles driven by
 * an autocatalytic reaction", Life 15 (2025), and their program growADE.c (Zenodo 15157714, CC-BY-4.0), re-implemented.
 *
 * - The volume follows osmosis at a finite rate: d(V/V0)/dt = lambda (A/A0) dC, with dC the excess of solutes and
 *   lipid inside over outside (every osmotic coefficient 1).
 * - The membrane is two leaflets of lipid area A+ and A-. Each exchanges lipid with the solution at rate rmem, towards
 *   the area its curvature asks for times the lipid in solution over C* (outside for the outer leaflet, inside for the
 *   inner one); lipid flips between them at rate rflip. The inner exchange changes the lipid in solution inside, Q
 *   being the lipid one area A0 holds per litre of the volume V0.
 * - Shape: a sphere when the volume fills the lipid area, else a capsule (two half-spheres of radius R joined by a
 *   cylinder of length L) of that area and volume. The area difference the shape asks for is 4 d M (M the integrated
 *   mean curvature).
 * - Elastic energy: curvature (kappa), area difference (beta kappa), stretching past the lipid area.
 * - Division, as the authors' energy rule: while the reduced volume v is below 1, at rate rdiv, a split into a sphere
 *   of a random share x < 1/2 of the volume and the rest is tried, the area difference shared in proportion to area;
 *   it happens when the two daughters' energy is below the parent's. One daughter is followed, the larger and the
 *   smaller alternately, as their program does. Or the equal rule: halving everything when v reaches 1/2.
 * - A sphere swollen past `burst_area` times its lipid area bursts, keeping `burst_kept` of its volume.
 * - A network's reactions run inside; each species may cross the membrane at mu (A/V) (outside - inside).
 *
 * Volumes and areas are relative to V0 and A0 (the sphere of radius R0), concentrations in mol/L, time in the units
 * of the rate constants. Integrated as the authors do, by classical Runge-Kutta of fourth order at a fixed step: the
 * capsule's radius is found by bisection to 1e-9, below which any adaptive error control would chase noise.
 */

typedef struct pl_capsule pl_capsule;

typedef struct pl_capsule_parameters {
    double lambda;           /* rate of the volume's change per unit area and excess concentration */
    double rmem, rflip;      /* lipid exchange with the solution, and flip-flop between leaflets */
    double rdiv;             /* rate of the attempts at division (energy rule) */
    int equal_division;      /* 1: halve everything when v reaches 1/2, instead of the energy rule */
    double cstar;            /* mol/L: lipid in solution in balance with a membrane */
    double elip;             /* mol/L: lipid in the outside solution */
    double R0, d, a;         /* m: reference radius, lipid length; m^2: lipid area */
    double kappa, beta, kappa_gauss;  /* bending modulus (kT units of the energy), area-difference factor, Gaussian */
    double burst_area;       /* A / A_lipid past which a swollen vesicle bursts */
    double burst_kept;       /* share of the volume kept at a burst */
    double largest;          /* V/V0 past which the run stops (the cell grew without bound) */
} pl_capsule_parameters;

/* The authors' standard set: lambda 1, rmem = rflip = 0.1, rdiv 10, C* 0.01 M, Elip 1.01 C*, R0 1 um, d 1 nm,
 * a 2e-19 m^2, kappa 1, beta 3/8, no Gaussian modulus, bursting at A = 1.1 A_lipid keeping 0.9 of the volume. */
pl_capsule_parameters pl_capsule_taneja_higgs(void);

typedef struct pl_capsule_shape {
    double volume, area;        /* V/V0, actual area / A0 */
    double lipid_area;          /* (A+ + A-) / 2 / A0 */
    double area_difference;     /* (A+ - A-) / A0 */
    double radius, length;      /* R/R0, L/R0 of the capsule (length 0 for a sphere) */
    double reduced;             /* v = V / V of the sphere of the lipid area */
    double s;                   /* lipid area over volume, relative to the sphere of radius R0 */
    double energy;              /* elastic energy, in the units of kappa */
    int swollen;                /* 1: a sphere stretched past its lipid area */
} pl_capsule_shape;

typedef enum pl_capsule_event_kind { PL_CAPSULE_DIVISION = 0, PL_CAPSULE_BURST } pl_capsule_event_kind;

typedef struct pl_capsule_event {
    pl_capsule_event_kind kind;
    double time;
    double reduced;   /* v just before */
    double volume;    /* V/V0 just before */
    double share;     /* division: the volume share of the daughter followed */
} pl_capsule_event;

/* A protocell of volume V0 and lipid area A0 at rest (leaflets as the sphere asks), with the network's initial
 * concentrations inside and C* of lipid. permeation and outside give, per network species, mu (0 for one the membrane
 * keeps) and the outside concentration (mol/L). The network must outlive the capsule; `seed` drives the attempts at
 * division. */
pl_status pl_capsule_create(const pl_capsule_parameters *parameters, const pl_network *network,
                            const double *permeation, const double *outside, uint64_t seed, pl_capsule **capsule,
                            pl_diagnostic *diagnostic);
void pl_capsule_free(pl_capsule *capsule);
/* The chemistry alone in a fixed compartment of A/V = A0/V0 (the authors' equilibration): over `duration`, step dt. */
pl_status pl_capsule_settle(pl_capsule *capsule, double duration, double dt, pl_diagnostic *diagnostic);
/* The protocell over `duration`, step dt: growth, exchange, chemistry, then attempts at division and bursts. It fails
 * with PL_ERROR_MODEL when the cell outgrows `largest` or becomes a tube thinner than R0 / 100 (the authors stop
 * there), leaving the state at that moment. */
pl_status pl_capsule_advance(pl_capsule *capsule, double duration, double dt, pl_diagnostic *diagnostic);

/* Sets one dynamic species' concentration inside (mol/L), as the authors' inactive cell, which starts without catalyst. */
pl_status pl_capsule_set_molar(pl_capsule *capsule, int species, double molar, pl_diagnostic *diagnostic);

double pl_capsule_time(const pl_capsule *capsule);
pl_capsule_shape pl_capsule_get_shape(const pl_capsule *capsule);
double pl_capsule_molar(const pl_capsule *capsule, int species);  /* inside, mol/L */
double pl_capsule_lipid(const pl_capsule *capsule);               /* lipid in solution inside, mol/L */
int pl_capsule_events(const pl_capsule *capsule, const pl_capsule_event **events);

#endif
