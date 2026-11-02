#ifndef PROTOLIFE_BALANCE_H
#define PROTOLIFE_BALANCE_H

/* A balance of matter a simulation checks at every step: what the system holds, less what came in, plus what left.
 * No model may change it; a step that does beyond the simulation's limit is undone and fails with
 * PL_ERROR_CONSERVATION. The changes below the limit add up in the drift. */
typedef struct pl_balance {
    const char *name;  /* "fatty acid", a solute's name, or an element's symbol ("C") */
    double value;      /* mol per litre of the vessel, as of the last step */
    double drift;      /* the sum of its changes over the steps so far, mol per litre of the vessel */
} pl_balance;

#endif
