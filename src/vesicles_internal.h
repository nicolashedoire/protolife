#ifndef PROTOLIFE_VESICLES_INTERNAL_H
#define PROTOLIFE_VESICLES_INTERNAL_H

#include "protolife/vesicles.h"

#define PL_PI 3.14159265358979323846
#define PL_AVOGADRO 6.02214076e23  /* 1/mol, exact since 2019 */
#define PL_BOLTZMANN 1.380649e-23  /* J/K, exact since 2019 */
#define PL_GAS 8.314462618         /* J/(mol K), N_A k */

typedef struct {
    char *name;
    double permeability, osmotic;
} solute_t;

typedef struct {
    double count, molecules, released;
    int generation, origin, inert, label;
    double coat;          /* mol/L of fed fatty acid on this record's membranes, entering fast */
    double intermediate;  /* mol/L of intermediate entering this record's membranes */
} record_t;

struct pl_vesicles {
    pl_vesicle_constants k;
    solute_t *solute;
    int solutes, solute_capacity, frozen;  /* frozen: solutes fixed once a vesicle exists */
    record_t *record;
    double *inside;                        /* max_records x solutes: mol of each solute in each vesicle */
    int records, max_records;
    double monomer;                        /* mol of free fatty acid per litre of the vessel */
    double *outside, *leaked, *feed, *renewed;
    double forming, newborn;               /* mol/L of intermediate making new vesicles, and fatty acid made into them */
    int nursery, feeds;
    double shear, divided, born, apart, washed, washout;
    long ruptures;
    double renewal, monomer_feed, renewed_monomer;
    int fatty_acid_held, lifeless;
    double anhydride, anhydride_feed, anhydride_fed;
    const pl_network *network;             /* reactions in every compartment, or NULL */
    int *species_of;                       /* per solute: the network species, or -1 */
    int *solute_of;                        /* per network species: the solute, or -1 */
    double *molar, *change;                /* per network species: working space of the reactions (mol/L) */
    double time, step, tolerance;
    pl_ode_statistics statistics;
    /* working space of one reading of the state */
    double *c_outside;                     /* mol/m^3 per solute */
    double *c_inside, *rate;               /* mol/m^3 per solute */
    pl_vesicle_shape *shape;               /* per record */
    double *work;                          /* 3 x state size: the Jacobian's shifted states, during a call */
    double fed;                            /* mol/L of fatty acid fed as micelles since the start */
    /* the balances checked at every step */
    double conservation;                   /* the largest relative change a step may make */
    pl_balance *balance;
    int balances, *balance_kind, *balance_index;  /* kind: 0 fatty acid, 1 a solute, 2 an element */
    double *before, *after, *scale, *scale_after;  /* per balance, during a call */
    /* the state before a step, to undo it */
    struct pl_vesicles *saved;
    record_t *saved_record;
    double *saved_inside, *saved_outside, *saved_leaked, *saved_renewed;
};

/* The layout of the integrated state, per litre of the vessel: each record's fatty acids (mol) and solutes inside
 * (mol); the free monomer and the solutes outside; each record's fatty acids released during the call (a ledger); the
 * coat and the intermediate of each record; the intermediate making vesicles and what it made; what the renewal
 * brought in, then what it took out (monomer, each solute); the anhydride. */
int pl_vesicles_state_size(const pl_vesicles *p);
void pl_vesicles_pack(const pl_vesicles *p, double *y);
/* The medium's volume per litre of the vessel and its osmolarity (osm/m^3), filling p->c_outside and p->shape. */
double pl_vesicles_read(pl_vesicles *p, const double *y, double *osmolarity);

/* The balances (vesicles_balances.c): made once the solutes are fixed; measured into value and scale; the state kept
 * before a step and put back. */
int pl_vesicles_prepare_balances(pl_vesicles *p);
void pl_vesicles_measure(const pl_vesicles *p, double *value, double *scale);
void pl_vesicles_save(pl_vesicles *p);
void pl_vesicles_restore(pl_vesicles *p);
void pl_vesicles_free_balances(pl_vesicles *p);

/* The end-of-step events (vesicles_events.c). */
void pl_vesicles_come_apart(pl_vesicles *p);
void pl_vesicles_births(pl_vesicles *p);
void pl_vesicles_divisions(pl_vesicles *p);
void pl_vesicles_wash(pl_vesicles *p, double dt);
int pl_vesicles_append(pl_vesicles *p, double count, double molecules, const double *inside, int label);

static inline double *pl_inside(const pl_vesicles *p, int record) { return p->inside + (size_t)record * (size_t)p->solutes; }

#endif
