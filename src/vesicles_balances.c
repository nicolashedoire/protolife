/* The balances a vesicle population checks at every step, and the saved state that undoes a step
 * (protolife/vesicles.h). */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "formula_internal.h"
#include "vesicles_internal.h"

enum { FATTY_ACID, SOLUTE, ELEMENT };

void pl_vesicles_free_balances(pl_vesicles *p) {
    free(p->balance), free(p->balance_kind), free(p->balance_index);
    free(p->before), free(p->after), free(p->scale), free(p->scale_after);
    free(p->saved), free(p->saved_record), free(p->saved_inside);
    free(p->saved_outside), free(p->saved_leaked), free(p->saved_renewed);
}

static int add_balance(pl_vesicles *p, const char *name, int kind, int index) {
    int b = p->balances++;
    p->balance[b] = (pl_balance){name, 0.0, 0.0};
    p->balance_kind[b] = kind, p->balance_index[b] = index;
    return b;
}

/* The elements some dynamic species holds and the network conserves. */
static void add_elements(pl_vesicles *p) {
    if (p->network == NULL) return;
    for (int e = 0; e < PL_ELEMENTS; e++) {
        const char *symbol = pl_element_symbol(e);
        int held = 0;
        for (int s = 0; s < p->solutes && !held; s++)
            held = p->species_of[s] >= 0 && pl_network_species_atoms(p->network, p->species_of[s], symbol) > 0;
        if (held && pl_network_conserves(p->network, symbol)) add_balance(p, symbol, ELEMENT, e);
    }
}

/* Once, before the first step: the solutes are then fixed, and so are the balances. 0, or -1 without memory. */
int pl_vesicles_prepare_balances(pl_vesicles *p) {
    if (p->saved != NULL) return 0;
    int solutes = p->solutes > 0 ? p->solutes : 1, most = 1 + solutes + PL_ELEMENTS;
    if (p->inside == NULL) {
        p->inside = calloc((size_t)p->max_records * (size_t)solutes, sizeof(double));
        if (p->inside == NULL) return -1;
    }
    p->frozen = 1;
    p->balance = calloc((size_t)most, sizeof *p->balance);
    p->balance_kind = calloc((size_t)most, sizeof *p->balance_kind);
    p->balance_index = calloc((size_t)most, sizeof *p->balance_index);
    p->before = calloc((size_t)most, sizeof *p->before), p->after = calloc((size_t)most, sizeof *p->after);
    p->scale = calloc((size_t)most, sizeof *p->scale), p->scale_after = calloc((size_t)most, sizeof *p->scale_after);
    p->saved = calloc(1, sizeof *p->saved);
    p->saved_record = calloc((size_t)p->max_records, sizeof *p->saved_record);
    p->saved_inside = calloc((size_t)p->max_records * (size_t)solutes, sizeof(double));
    double **arrays[] = {&p->saved_outside, &p->saved_leaked, &p->saved_renewed};
    for (size_t i = 0; i < sizeof arrays / sizeof *arrays; i++) *arrays[i] = calloc((size_t)solutes, sizeof(double));
    if (p->balance == NULL || p->balance_kind == NULL || p->balance_index == NULL || p->before == NULL ||
        p->after == NULL || p->scale == NULL || p->scale_after == NULL || p->saved == NULL || p->saved_record == NULL || p->saved_inside == NULL ||
        p->saved_outside == NULL || p->saved_leaked == NULL || p->saved_renewed == NULL) {
        pl_vesicles_free_balances(p);
        p->balance = NULL, p->balance_kind = p->balance_index = NULL, p->before = p->after = p->scale = p->scale_after = NULL;
        p->saved = NULL, p->saved_record = NULL;
        p->saved_inside = p->saved_outside = p->saved_leaked = p->saved_renewed = NULL;
        p->balances = 0;
        return -1;
    }
    add_balance(p, "fatty acid", FATTY_ACID, 0);
    for (int s = 0; s < p->solutes; s++)
        if (p->species_of[s] < 0) add_balance(p, p->solute[s].name, SOLUTE, s);
    add_elements(p);
    pl_vesicles_measure(p, p->before, p->scale);
    for (int b = 0; b < p->balances; b++) p->balance[b].value = p->before[b];
    return 0;
}

/* A solute's balance: what the vessel holds less what the renewal brought in (net); the size of its terms. */
static double solute_balance(const pl_vesicles *p, int s, double *scale) {
    double total = pl_vesicles_solute_total(p, s);
    *scale = fabs(total) + fabs(p->renewed[s]);
    return total - p->renewed[s];
}

void pl_vesicles_measure(const pl_vesicles *p, double *value, double *scale) {
    for (int b = 0; b < p->balances; b++) {
        int kind = p->balance_kind[b], index = p->balance_index[b];
        if (kind == FATTY_ACID) {
            double in_vessel = pl_vesicles_fatty_acid(p) + 2.0 * p->anhydride;
            value[b] = in_vessel + p->washed - p->renewed_monomer - p->fed - 2.0 * p->anhydride_fed;
            scale[b] = in_vessel + p->washed + fabs(p->renewed_monomer) + p->fed + 2.0 * p->anhydride_fed;
        } else if (kind == SOLUTE) {
            value[b] = solute_balance(p, index, &scale[b]);
        } else {
            const char *symbol = pl_element_symbol(index);
            value[b] = scale[b] = 0.0;
            for (int s = 0; s < p->solutes; s++) {
                int k = p->species_of[s], atoms = k >= 0 ? pl_network_species_atoms(p->network, k, symbol) : 0;
                if (atoms == 0) continue;
                double size;
                value[b] += atoms * solute_balance(p, s, &size);
                scale[b] += atoms * size;
            }
        }
    }
}

void pl_vesicles_save(pl_vesicles *p) {
    *p->saved = *p;
    size_t solutes = (size_t)p->solutes;
    memcpy(p->saved_record, p->record, (size_t)p->records * sizeof *p->record);
    memcpy(p->saved_inside, p->inside, (size_t)p->records * solutes * sizeof(double));
    memcpy(p->saved_outside, p->outside, solutes * sizeof(double));
    memcpy(p->saved_leaked, p->leaked, solutes * sizeof(double));
    memcpy(p->saved_renewed, p->renewed, solutes * sizeof(double));
}

/* Back to the saved state; the work done in the step undone still counts in the statistics. */
void pl_vesicles_restore(pl_vesicles *p) {
    pl_ode_statistics work = p->statistics;
    *p = *p->saved;
    size_t solutes = (size_t)p->solutes;
    memcpy(p->record, p->saved_record, (size_t)p->records * sizeof *p->record);
    memcpy(p->inside, p->saved_inside, (size_t)p->records * solutes * sizeof(double));
    memcpy(p->outside, p->saved_outside, solutes * sizeof(double));
    memcpy(p->leaked, p->saved_leaked, solutes * sizeof(double));
    memcpy(p->renewed, p->saved_renewed, solutes * sizeof(double));
    p->statistics = work;
}
