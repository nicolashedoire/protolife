#ifndef PROTOLIFE_VESICLES_H
#define PROTOLIFE_VESICLES_H

#include "protolife/balance.h"
#include "protolife/diagnostic.h"
#include "protolife/network.h"
#include "protolife/ode.h"
#include "protolife/rodas.h"

/*
 * A population of fatty acid vesicles in a stirred medium: each one a membrane of fatty acids around its own water and
 * solutes. A record stands for many identical vesicles (per litre of the vessel). After Genesis world/vesicles.c.
 *
 * - Water follows osmosis at once: a vesicle takes the volume that matches the osmolarities inside and out, ideal
 *   (van 't Hoff), no larger than the sphere its membrane can close. A vesicle that would hold more is a sphere under
 *   tension, Laplace's 2 tau / r = delta Pi.
 * - Fatty acids leave a membrane at the measured desorption rate and come back in proportion to the free monomer; a
 *   relaxed membrane is in balance with the critical vesicle concentration. Tension holds them: leaving is slowed by
 *   exp(-tau a / 2kT), tau capped at the rupture tension, so a tense vesicle grows at the expense of relaxed ones
 *   (Chen, Roberts and Szostak, Science 305, 1474, 2004).
 * - Each solute crosses at its permeability, J = P A (c_out - c_in) (Sacerdote and Szostak, PNAS 102, 6004, 2005).
 * - Past the rupture tension the membrane opens, lets its contents out in the same proportion until the tension is back
 *   at the rupture value, and reseals (Chen 2004b).
 * - Micelles fed to the medium (Chen and Szostak, Biophys. J. 87, 988, 2004): up to 42 % of the membranes' fatty acid
 *   coats them and enters fast; the rest forms an intermediate that enters slowly, shared between the membranes and new
 *   vesicles of 45 nm.
 * - With a network attached, its reactions run in every vesicle and in the medium, each at its own concentrations; the
 *   medium may be renewed while the vesicles are held back, as behind a dialysis membrane.
 * - An optional oil layer of oleic anhydride (Walde et al., J. Am. Chem. Soc. 116, 11649, 1994) hydrolyses at its face
 *   and in the membranes.
 * - A vesicle shrinking below `apart_radius` comes apart; the vessel may be washed through (a stream).
 * - Under a shear at least the critical one (Zhu and Szostak, J. Am. Chem. Soc. 131, 5705, 2009), a vesicle stretched
 *   thin (reduced volume below a threshold) divides into the fewest equal daughters whose reduced volume reaches it.
 *
 * Each call of pl_vesicles_advance integrates the continuous part with RODAS3 (protolife/rodas.h) and Genesis's
 * Jacobian by differences, as Genesis does: the switch between a relaxed vesicle and a tense sphere is a kink where a
 * Newton iteration stalls. Ruptures, comings apart, births, divisions and the washout are then applied, in that order.
 * A call that fails changes nothing.
 *
 * Every call also checks the balances that must hold whatever the model does: the fatty acid (in membranes, free,
 * entering membranes, as anhydride; with what was fed as micelles or poured as anhydride, renewed and washed out), each
 * solute outside the network (with what was renewed and washed out), and each element the network conserves
 * (pl_network_conserves: carbon for a sugar network, not hydrogen or oxygen, which water supplies). A step that changes
 * one of them by more than `conservation` times the size of its terms is undone, and the call fails with
 * PL_ERROR_CONSERVATION naming the balance. The changes below that limit add up in each balance's drift.
 * Amounts are per litre of the vessel, concentrations in the medium in mol/L, lengths in m, time in s.
 */

typedef struct pl_vesicles pl_vesicles;

typedef struct pl_vesicle_constants {
    double area_per_molecule;     /* m^2, one fatty acid in one leaflet */
    double desorption;            /* 1/s, from a relaxed membrane */
    double cvc;                   /* mol/L, free fatty acid in balance with relaxed membranes */
    double rupture;               /* N/m */
    double kelvin;
    double coat_rate, slow_rate;  /* 1/s: the fast and slow entries of fed micelles */
    double coat_share;            /* the coat a membrane takes at once, a share of its fatty acid */
    double nucleation;            /* kappa: how the intermediate shares itself with new vesicles */
    double born_radius;           /* m, of the vesicles the intermediate makes */
    double critical_shear;        /* 1/s */
    double division_reduced;      /* the reduced volume below which a sheared vesicle divides */
    double least_radius;          /* m: no daughter smaller */
    double apart_radius;          /* m: a vesicle shrinking below it comes apart */
    double interface_hydrolysis;  /* mol/L/s: anhydride hydrolysed at its layer's face (0: none) */
    double membrane_hydrolysis;   /* 1/s: of the anhydride the membranes hold */
    double uptake;                /* anhydride a membrane holds at most, per fatty acid in it (0: none) */
    double closing;               /* 1/s: monomer above the cvc closing into new vesicles (0: none) */
} pl_vesicle_constants;

/* Oleic acid / oleate at pH 8.5 in 0.2 M bicine, about 23 °C: area 0.32 nm^2, desorption 3.0/min (Chen 2004b), cvc
 * 0.5 mM, rupture 10 dyn/cm; micelle entry 2.5/s and 0.1/s, coat 42 %, new vesicles of 45 nm (Chen and Szostak 2004);
 * critical shear 15/s (Zhu and Szostak 2009); division below a reduced volume of 0.55 (assumed); daughters of 20 nm at
 * least and coming apart below 10 nm (assumed). The anhydride and the closing are off. */
pl_vesicle_constants pl_vesicle_oleate(void);

typedef struct pl_vesicle_shape {
    double volume, area, radius, tension;  /* m^3, m^2, m (of the area's sphere), N/m */
    int spherical;                         /* 1 when the membrane is stretched round its contents */
} pl_vesicle_shape;

/* One record: `count` identical vesicles per litre of the vessel. */
typedef struct pl_vesicle {
    double count;           /* vesicles per litre of the vessel */
    double molecules;       /* fatty acids in each membrane */
    double released;        /* fatty acids that left each membrane since the start */
    int generation;         /* divisions since the record's vesicles were made */
    int origin;             /* 0 for the vesicles of the start, k for those born of the k-th feed of micelles */
    int inert;              /* 1: no reaction inside (a control of the same size and contents) */
    int label;              /* the caller's name for the record */
    double coat;            /* mol/L of fed fatty acid on the record's membranes, entering fast */
    double intermediate;    /* mol/L of intermediate entering the record's membranes, slowly */
    const double *inside;   /* mol of each solute in each vesicle; valid until the population changes */
} pl_vesicle;

/* Counts since the start, per litre of the vessel. */
typedef struct pl_vesicle_ledger {
    long ruptures;
    int feeds;               /* feeds of micelles */
    double divided, born;    /* vesicles made by divisions and by micelles */
    double apart;            /* vesicles that came apart */
    double washed;           /* mol of fatty acid that left in vesicles (washout, dilution) */
    double renewed_monomer;  /* mol of free fatty acid brought in, net, by the renewal and dilutions */
    double fed;              /* mol of fatty acid fed as micelles */
    const double *renewed;   /* mol of each solute brought in, net (what left inside vesicles counted as leaving) */
    const double *leaked;    /* mol of each solute let out by ruptures */
} pl_vesicle_ledger;

/* An empty population: no solute, no vesicle; at most `max_records` records (0 for 24). */
pl_status pl_vesicles_create(const pl_vesicle_constants *constants, int max_records, pl_vesicles **population,
                             pl_diagnostic *diagnostic);
void pl_vesicles_free(pl_vesicles *population);

/* A solute of the medium (before any vesicle): its permeability (m/s, 0 for one the membrane keeps) and osmoles per
 * mole. Its index goes to *index when index is not NULL. */
pl_status pl_vesicles_add_solute(pl_vesicles *population, const char *name, double permeability, double osmotic,
                                 int *index, pl_diagnostic *diagnostic);
/* A network's reactions in every vesicle and in the medium (before any vesicle). Its dynamic species become solutes
 * (osmotic coefficient 1), crossing at permeability(network, species, context) m/s; reservoir and implicit species keep
 * their fixed concentrations everywhere. The network must outlive the population. */
pl_status pl_vesicles_attach_network(pl_vesicles *population, const pl_network *network,
                                     double (*permeability)(const pl_network *network, int species, void *context),
                                     void *context, pl_diagnostic *diagnostic);
/* The solute standing for a network species, or -1; the network species a solute stands for, or -1. */
int pl_vesicles_solute_of(const pl_vesicles *population, int species);
int pl_vesicles_species_of(const pl_vesicles *population, int solute);

/* A record of `count` vesicles per litre, each of `molecules` fatty acids holding `inside` (mol of each solute, NULL
 * for none), under the caller's `label`. Its index goes to *record when record is not NULL. */
pl_status pl_vesicles_add(pl_vesicles *population, double count, double molecules, const double *inside, int label,
                          int *record, pl_diagnostic *diagnostic);
pl_status pl_vesicles_set_inert(pl_vesicles *population, int record, int inert, pl_diagnostic *diagnostic);

/* The medium. Amounts are mol per litre of the vessel; feeds are mol/L of the fresh medium. */
pl_status pl_vesicles_set_outside(pl_vesicles *population, int solute, double amount, pl_diagnostic *diagnostic);
pl_status pl_vesicles_set_monomer(pl_vesicles *population, double amount, pl_diagnostic *diagnostic);
pl_status pl_vesicles_set_feed(pl_vesicles *population, int solute, double molar, pl_diagnostic *diagnostic);
pl_status pl_vesicles_set_monomer_feed(pl_vesicles *population, double molar, pl_diagnostic *diagnostic);
/* The share of the medium replaced each second (1/s), the vesicles held back; with fatty_acid_held, the fatty acid
 * is held back too, as in a dialysis bag. */
pl_status pl_vesicles_set_renewal(pl_vesicles *population, double rate, int fatty_acid_held, pl_diagnostic *diagnostic);
pl_status pl_vesicles_set_shear(pl_vesicles *population, double rate, pl_diagnostic *diagnostic);
pl_status pl_vesicles_set_washout(pl_vesicles *population, double rate, pl_diagnostic *diagnostic);
/* 1: no reaction in any vesicle, newborns included (a control vessel); the medium still reacts. */
void pl_vesicles_set_lifeless(pl_vesicles *population, int lifeless);
/* The oil layer: mol of anhydride per litre of the vessel, and what is poured on it (mol/L/s). */
pl_status pl_vesicles_set_anhydride(pl_vesicles *population, double amount, double feed, pl_diagnostic *diagnostic);
/* The relative tolerance of the integration (default 1e-8). */
pl_status pl_vesicles_set_tolerance(pl_vesicles *population, double relative, pl_diagnostic *diagnostic);

/* Over dt (s). Ruptures, comings apart, births, divisions and the washout are applied at the end of the step. */
pl_status pl_vesicles_advance(pl_vesicles *population, double dt, pl_diagnostic *diagnostic);
/* Micelles of `amount` mol/L of fatty acid fed to the medium: coat, intermediate and new vesicles. */
void pl_vesicles_feed_micelles(pl_vesicles *population, double amount);
/* A share of everything but the oil layer replaced at once by fresh medium: vesicles and micelles leave (counted as
 * washed), the medium's solutes and monomer are renewed. */
void pl_vesicles_dilute(pl_vesicles *population, double share);
/* Ruptures now, for a change of the medium made between steps. */
void pl_vesicles_rupture(pl_vesicles *population);

/* The balances checked, as of the last call of pl_vesicles_advance (none before the first): their number, and the list
 * in *balances, valid until the population changes. */
int pl_vesicles_balances(const pl_vesicles *population, const pl_balance **balances);
/* The largest change of a balance one step may make, relative to the size of its terms (default 1e-9). */
pl_status pl_vesicles_set_conservation(pl_vesicles *population, double relative, pl_diagnostic *diagnostic);

double pl_vesicles_time(const pl_vesicles *population);
int pl_vesicles_records(const pl_vesicles *population);
int pl_vesicles_solutes(const pl_vesicles *population);
const char *pl_vesicles_solute_name(const pl_vesicles *population, int solute);
pl_vesicle pl_vesicles_record(const pl_vesicles *population, int record);
/* The record of this label, or -1 (none, or it came apart). */
int pl_vesicles_find(const pl_vesicles *population, int label);
pl_vesicle_shape pl_vesicles_shape(const pl_vesicles *population, int record);
/* The medium's volume per litre of the vessel (what the vesicles do not hold) and its osmolarity (osm/m^3). */
double pl_vesicles_medium(const pl_vesicles *population, double *osmolarity);
double pl_vesicles_outside(const pl_vesicles *population, int solute);  /* mol per litre of the vessel */
double pl_vesicles_monomer(const pl_vesicles *population);              /* mol per litre of the vessel */
double pl_vesicles_shear(const pl_vesicles *population);
/* mol/L of fatty acid on its way into the membranes: the coat, the intermediate entering them, the one making new
 * vesicles and what it made, summed over the records. */
double pl_vesicles_entering(const pl_vesicles *population);
/* mol/L of the intermediate making new vesicles, and of the fatty acid it made waiting to close into them. */
double pl_vesicles_forming(const pl_vesicles *population);
double pl_vesicles_membranes(const pl_vesicles *population);           /* mol/L of fatty acid in membranes */
double pl_vesicles_fatty_acid(const pl_vesicles *population);          /* all fatty acid, mol/L */
double pl_vesicles_solute_total(const pl_vesicles *population, int solute);
pl_vesicle_ledger pl_vesicles_ledger(const pl_vesicles *population);
void pl_vesicles_get_statistics(const pl_vesicles *population, pl_ode_statistics *statistics);

/* The fatty acids a sphere of radius r needs; the shape of a membrane of `molecules` holding `osmoles` in a medium of
 * `osmolarity` (osm/m^3). */
double pl_vesicle_molecules_for_radius(const pl_vesicle_constants *constants, double radius);
pl_vesicle_shape pl_vesicle_shape_of(const pl_vesicle_constants *constants, double molecules, double osmoles,
                                     double osmolarity);

#endif
