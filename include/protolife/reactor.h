#ifndef PROTOLIFE_REACTOR_H
#define PROTOLIFE_REACTOR_H

#include "protolife/diagnostic.h"
#include "protolife/network.h"
#include "protolife/ode.h"

/*
 * A well-stirred vessel running a network's rate equations, d[S]/dt = sum over reactions of (net coefficient) * rate,
 * integrated by pl_ode with concentrations kept nonnegative. Closed, or open as a continuously stirred tank reactor
 * when `dilution` > 0: every dynamic species then also flows out at dilution * [S] and comes in at
 * dilution * feed[S]. Reservoir and implicit species stay exactly at their fixed concentrations. Concentrations in
 * mol/L, time in s.
 */

typedef struct pl_reactor pl_reactor;

typedef struct pl_reactor_options {
    double relative_tolerance;  /* 0 for 1e-8 */
    double absolute_tolerance;  /* mol/L, the size under which a concentration is noise; 0 for 1e-14 */
    double dilution;            /* 1/s, >= 0; 0 for a closed vessel */
    const double *feed;         /* mol/L for every species (only dynamic ones are used); NULL for none */
    long max_steps;             /* per call of pl_reactor_advance; 0 for 100000 */
} pl_reactor_options;

/* A reactor at time 0 with the network's initial and fixed concentrations. The network must outlive it and must not
 * change while it runs. options may be NULL for a closed vessel with the default tolerances. */
pl_status pl_reactor_create(const pl_network *network, const pl_reactor_options *options, pl_reactor **reactor,
                            pl_diagnostic *diagnostic);
/* Integrates up to t_end (s); on failure the reactor keeps its time and concentrations (protolife/ode.h). */
pl_status pl_reactor_advance(pl_reactor *reactor, double t_end, pl_diagnostic *diagnostic);
/* Sets one species' concentration (dynamic, or the fixed value of a reservoir) and restarts the integration. */
pl_status pl_reactor_set_molar(pl_reactor *reactor, int species, double molar, pl_diagnostic *diagnostic);
double pl_reactor_time(const pl_reactor *reactor);
/* The concentration of every species (mol/L), owned by the reactor and valid until it is freed. */
const double *pl_reactor_molar(const pl_reactor *reactor);
void pl_reactor_get_statistics(const pl_reactor *reactor, pl_ode_statistics *statistics);
void pl_reactor_free(pl_reactor *reactor);

#endif
