#ifndef PROTOLIFE_NETWORK_H
#define PROTOLIFE_NETWORK_H

#include "protolife/diagnostic.h"

/*
 * A reaction network: species with a formula, a charge and a role, and reactions with a mass-action rate law,
 * r = k * product over the left-hand side of [X]^coefficient, in mol/L/s with concentrations in mol/L (the unit of
 * published rate constants).
 *
 * Roles: dynamic species change; reservoir species are held at a fixed concentration from outside (an open system);
 * implicit species (H+ at a fixed pH, water) are folded into the rate constants and count only in the balance of atoms
 * and charge. A reaction may carry implicit terms, which count in the balance only; a lumped reaction stands for
 * several steps and is exempt from the balance check. Every reaction cites its source.
 *
 * The text format (Genesis's .net files), one statement per line, '#' starting a comment:
 *   temperature <kelvin>                      (298.15 when absent)
 *   species <name> <formula> <charge> dynamic|reservoir|implicit <mol/L>
 *   reaction <equation> ; k=<value> [; implicit=<equation>] [; lumped] ; source=<text>
 * An equation is "2 X + Y -> Z"; a side may be empty; a coefficient is a positive integer before a name.
 */

typedef struct pl_network pl_network;

typedef enum pl_role { PL_DYNAMIC = 0, PL_RESERVOIR, PL_IMPLICIT } pl_role;

/* An empty network at 298.15 K. */
pl_status pl_network_create(pl_network **network, pl_diagnostic *diagnostic);
/* A network from text in the format above. On failure the diagnostic's code is the line number. */
pl_status pl_network_parse(const char *text, pl_network **network, pl_diagnostic *diagnostic);
/* A network from a file in the format above. */
pl_status pl_network_load(const char *path, pl_network **network, pl_diagnostic *diagnostic);
void pl_network_free(pl_network *network);

/* A species: its name (no blank, not starting with a digit, unique), its formula ("CH2(COOH)2"), its charge, its role
 * and its initial (dynamic) or fixed (reservoir, implicit) concentration in mol/L. */
pl_status pl_network_add_species(pl_network *network, const char *name, const char *formula, int charge, pl_role role,
                                 double molar, pl_diagnostic *diagnostic);
/* A reaction: its equation over species already added, its rate constant (>= 0, in mol/L and s), its implicit terms
 * as an equation (NULL for none), whether it is lumped, and its source (required). */
pl_status pl_network_add_reaction(pl_network *network, const char *equation, double rate, const char *implicit,
                                  int lumped, const char *source, pl_diagnostic *diagnostic);
pl_status pl_network_set_temperature(pl_network *network, double kelvin, pl_diagnostic *diagnostic);

/* Every reaction that is not lumped conserves each element and the charge, implicit terms included; otherwise
 * PL_ERROR_CONSERVATION, naming the first reaction and the first quantity that differ. */
pl_status pl_network_check_balance(const pl_network *network, pl_diagnostic *diagnostic);

int pl_network_species_count(const pl_network *network);
int pl_network_reaction_count(const pl_network *network);
double pl_network_temperature(const pl_network *network);
/* The index of a species, or -1. */
int pl_network_find(const pl_network *network, const char *name);
const char *pl_network_species_name(const pl_network *network, int species);
const char *pl_network_species_formula(const pl_network *network, int species);
int pl_network_species_charge(const pl_network *network, int species);
pl_role pl_network_species_role(const pl_network *network, int species);
double pl_network_species_molar(const pl_network *network, int species);
double pl_network_reaction_constant(const pl_network *network, int reaction);
const char *pl_network_reaction_source(const pl_network *network, int reaction);
int pl_network_reaction_lumped(const pl_network *network, int reaction);

/* The rate of every reaction (mol/L/s) at the concentrations of every species (mol/L). */
void pl_network_rates(const pl_network *network, const double *molar, double *rate);
/* d[S]/dt of every species from the reactions alone (mol/L/s); 0 for reservoir and implicit species. */
void pl_network_derivative(const pl_network *network, const double *molar, double *dmolar);
/* The Jacobian of pl_network_derivative, row-major species x species (1/s); columns of fixed species are 0. */
void pl_network_jacobian(const pl_network *network, const double *molar, double *jacobian);

#endif
