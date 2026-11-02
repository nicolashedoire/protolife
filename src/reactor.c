/* A well-stirred reactor over pl_ode (protolife/reactor.h), after Genesis's chemistry/mean_field.c. */
#include "protolife/reactor.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "diagnostic_internal.h"
#include "formula_internal.h"

struct pl_reactor {
    const pl_network *network;
    pl_ode *ode;
    int n, dynamic;     /* species, and those integrated */
    int *index;         /* the species of each integrated component */
    double dilution;
    double *feed;       /* every species */
    double *molar;      /* every species: the state the caller sees */
    double *scratch;    /* every species: concentrations inside a callback */
    double *derivative; /* every species */
    double *jacobian;   /* species x species */
    /* one balance per conserved element: two ledgers each after the species in the state (atoms in, atoms out) */
    int elements, *element, *atoms;  /* atoms: elements x dynamic, per integrated component */
    pl_balance *balance;
    double conservation, *before, *scale, *saved;
};

/* Every species' concentration from the integrated components y, the fixed ones from `molar`. */
static void expand(pl_reactor *r, const double *y) {
    memcpy(r->scratch, r->molar, (size_t)r->n * sizeof *r->scratch);
    for (int i = 0; i < r->dynamic; i++) r->scratch[r->index[i]] = y[i];
}

static int rhs(double t, const double *y, double *dydt, void *context) {
    (void)t;
    pl_reactor *r = context;
    expand(r, y);
    pl_network_derivative(r->network, r->scratch, r->derivative);
    for (int i = 0; i < r->dynamic; i++) {
        int s = r->index[i];
        dydt[i] = r->derivative[s] + r->dilution * (r->feed[s] - y[i]);
    }
    for (int e = 0; e < r->elements; e++) {  /* the atoms the feed brings in and the outflow takes away */
        double in = 0.0, out = 0.0;
        for (int i = 0; i < r->dynamic; i++) {
            int a = r->atoms[e * r->dynamic + i];
            in += a * r->dilution * r->feed[r->index[i]], out += a * r->dilution * y[i];
        }
        dydt[r->dynamic + 2 * e] = in, dydt[r->dynamic + 2 * e + 1] = out;
    }
    return 0;
}

static int jacobian(double t, const double *y, double *out, void *context) {
    (void)t;
    pl_reactor *r = context;
    int size = r->dynamic + 2 * r->elements;
    expand(r, y);
    pl_network_jacobian(r->network, r->scratch, r->jacobian);
    memset(out, 0, (size_t)size * (size_t)size * sizeof *out);
    for (int i = 0; i < r->dynamic; i++)
        for (int j = 0; j < r->dynamic; j++)
            out[i * size + j] = r->jacobian[r->index[i] * r->n + r->index[j]] - (i == j ? r->dilution : 0.0);
    for (int e = 0; e < r->elements; e++)
        for (int j = 0; j < r->dynamic; j++) out[(r->dynamic + 2 * e + 1) * size + j] = r->dilution * r->atoms[e * r->dynamic + j];
    return 0;
}

/* Each balance: the atoms held, less those brought in, plus those taken away; and the size of its terms. */
static void measure(const pl_reactor *r, const double *y, double *value, double *scale) {
    for (int e = 0; e < r->elements; e++) {
        double held = 0.0, in = y[r->dynamic + 2 * e], out = y[r->dynamic + 2 * e + 1];
        for (int i = 0; i < r->dynamic; i++) held += r->atoms[e * r->dynamic + i] * y[i];
        value[e] = held - in + out, scale[e] = fabs(held) + fabs(in) + fabs(out);
    }
}

/* The elements some dynamic species holds and the network conserves, and their atoms per component. */
static int find_elements(pl_reactor *r) {
    int *element = malloc(PL_ELEMENTS * sizeof *element), count = 0;
    if (element == NULL) return -1;
    for (int e = 0; e < PL_ELEMENTS; e++) {
        const char *symbol = pl_element_symbol(e);
        int held = 0;
        for (int i = 0; i < r->dynamic && !held; i++) held = pl_network_species_atoms(r->network, r->index[i], symbol) > 0;
        if (held && pl_network_conserves(r->network, symbol)) element[count++] = e;
    }
    r->element = element, r->elements = count;
    int slots = count > 0 ? count : 1;
    r->atoms = malloc((size_t)slots * (size_t)r->dynamic * sizeof *r->atoms);
    r->balance = calloc((size_t)slots, sizeof *r->balance);
    r->before = malloc(4 * (size_t)slots * sizeof *r->before);
    if (r->atoms == NULL || r->balance == NULL || r->before == NULL) return -1;
    r->scale = r->before + slots;
    for (int e = 0; e < count; e++) {
        r->balance[e].name = pl_element_symbol(element[e]);
        for (int i = 0; i < r->dynamic; i++)
            r->atoms[e * r->dynamic + i] = pl_network_species_atoms(r->network, r->index[i], r->balance[e].name);
    }
    return 0;
}

void pl_reactor_free(pl_reactor *r) {
    if (r == NULL) return;
    pl_ode_free(r->ode);
    free(r->index), free(r->feed), free(r->molar), free(r->scratch), free(r->derivative), free(r->jacobian);
    free(r->element), free(r->atoms), free(r->balance), free(r->before), free(r->saved);
    free(r);
}

static pl_status check(const pl_network *network, const pl_reactor_options *o, pl_diagnostic *diagnostic) {
    if (network == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the network is required");
    if (pl_network_species_count(network) == 0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the network has no species");
    if (o == NULL) return PL_OK;
    if (!(o->relative_tolerance >= 0.0) || !(o->absolute_tolerance >= 0.0) || !isfinite(o->relative_tolerance) ||
        !isfinite(o->absolute_tolerance) || o->max_steps < 0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "tolerances and max_steps must be finite and >= 0");
    if (!(o->dilution >= 0.0) || !isfinite(o->dilution))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the dilution must be finite and >= 0");
    if (o->feed != NULL)
        for (int s = 0; s < pl_network_species_count(network); s++)
            if (!(o->feed[s] >= 0.0) || !isfinite(o->feed[s]))
                return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the feed of %s must be finite and >= 0",
                               pl_network_species_name(network, s));
    return PL_OK;
}

pl_status pl_reactor_create(const pl_network *network, const pl_reactor_options *options, pl_reactor **out,
                            pl_diagnostic *diagnostic) {
    if (out == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the reactor's address is required");
    *out = NULL;
    pl_status status = check(network, options, diagnostic);
    if (status != PL_OK) return status;
    pl_reactor_options o = options != NULL ? *options : (pl_reactor_options){0};
    pl_reactor *r = calloc(1, sizeof *r);
    if (r == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, 0.0, "the reactor could not be allocated");
    int n = pl_network_species_count(network);
    size_t bytes = (size_t)n * sizeof(double);
    r->network = network, r->n = n, r->dilution = o.dilution;
    r->index = malloc((size_t)n * sizeof *r->index);
    r->feed = calloc((size_t)n, sizeof *r->feed);
    r->molar = malloc(bytes), r->scratch = malloc(bytes), r->derivative = malloc(bytes);
    r->jacobian = malloc(bytes * (size_t)n);
    double *y0 = malloc(bytes), *absolute = malloc(bytes);
    if (r->index == NULL || r->feed == NULL || r->molar == NULL || r->scratch == NULL || r->derivative == NULL ||
        r->jacobian == NULL || y0 == NULL || absolute == NULL) {
        free(y0), free(absolute), pl_reactor_free(r);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, 0.0, "the reactor could not be allocated");
    }
    for (int s = 0; s < n; s++) {
        r->molar[s] = pl_network_species_molar(network, s);
        if (o.feed != NULL) r->feed[s] = o.feed[s];
        if (pl_network_species_role(network, s) == PL_DYNAMIC) r->index[r->dynamic++] = s;
    }
    if (r->dynamic == 0) {
        free(y0), free(absolute), pl_reactor_free(r);
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, 0.0, "the network has no dynamic species to integrate");
    }
    r->conservation = 1e-9;
    double floor = o.absolute_tolerance > 0.0 ? o.absolute_tolerance : 1e-14;
    int size = r->dynamic;
    if (find_elements(r) == 0) {
        size = r->dynamic + 2 * r->elements;
        double *bigger_y = realloc(y0, (size_t)size * sizeof *y0), *bigger_a = bigger_y != NULL ? realloc(absolute, (size_t)size * sizeof *absolute) : NULL;
        if (bigger_y != NULL) y0 = bigger_y;
        if (bigger_a != NULL) absolute = bigger_a;
        r->saved = malloc((size_t)size * sizeof *r->saved);
        if (bigger_y == NULL || bigger_a == NULL || r->saved == NULL) size = -1;
    } else {
        size = -1;
    }
    if (size < 0) {
        free(y0), free(absolute), pl_reactor_free(r);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, 0.0, "the reactor's balances could not be allocated");
    }
    for (int i = 0; i < r->dynamic; i++) y0[i] = r->molar[r->index[i]], absolute[i] = floor;
    for (int e = 0; e < r->elements; e++) {  /* ledgers start at zero, kept to the accuracy of the atoms they count */
        double atoms = 0.0;
        for (int i = 0; i < r->dynamic; i++) atoms += r->atoms[e * r->dynamic + i];
        y0[r->dynamic + 2 * e] = y0[r->dynamic + 2 * e + 1] = 0.0;
        absolute[r->dynamic + 2 * e] = absolute[r->dynamic + 2 * e + 1] = floor * (atoms > 0.0 ? atoms : 1.0);
    }
    measure(r, y0, r->before, r->scale);
    for (int e = 0; e < r->elements; e++) r->balance[e].value = r->before[e];
    {
        pl_ode_options ode = {.size = size, .rhs = rhs, .jacobian = jacobian, .context = r,
                              .relative_tolerance = o.relative_tolerance > 0.0 ? o.relative_tolerance : 1e-8,
                              .absolute_tolerance = absolute, .nonnegative = 1,
                              .max_steps = o.max_steps > 0 ? o.max_steps : 100000};
        status = pl_ode_create(&ode, 0.0, y0, &r->ode, diagnostic);
    }
    free(y0), free(absolute);
    if (status != PL_OK) {
        pl_reactor_free(r);
        return status;
    }
    *out = r;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_reactor_advance(pl_reactor *r, double t_end, pl_diagnostic *diagnostic) {
    if (r == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the reactor is required");
    int size = pl_ode_size(r->ode);
    double t0 = pl_ode_time(r->ode), *after = r->before + 2 * (r->elements > 0 ? r->elements : 1), *scale_after = after + (r->elements > 0 ? r->elements : 1);
    memcpy(r->saved, pl_ode_state(r->ode), (size_t)size * sizeof *r->saved);
    measure(r, r->saved, r->before, r->scale);
    pl_status status = pl_ode_advance(r->ode, t_end, NULL, diagnostic);
    if (status != PL_OK) return status;
    const double *y = pl_ode_state(r->ode);
    measure(r, y, after, scale_after);
    for (int e = 0; e < r->elements; e++) {
        double change = after[e] - r->before[e], limit = r->conservation * fmax(r->scale[e], scale_after[e]);
        if (fabs(change) <= limit) continue;
        pl_ode_restart(r->ode, t0, r->saved, NULL);
        return pl_fail(diagnostic, PL_ERROR_CONSERVATION, e, t_end,
                       "the balance of %s changed by %.3g mol/L (%.3g allowed) between t = %.17g and %.17g s: the step is "
                       "undone", r->balance[e].name, change, limit, t0, t_end);
    }
    for (int e = 0; e < r->elements; e++) r->balance[e].drift += after[e] - r->before[e], r->balance[e].value = after[e];
    for (int i = 0; i < r->dynamic; i++) r->molar[r->index[i]] = y[i];
    return PL_OK;
}

int pl_reactor_balances(const pl_reactor *r, const pl_balance **balances) {
    if (balances != NULL) *balances = r != NULL ? r->balance : NULL;
    return r != NULL ? r->elements : 0;
}

pl_status pl_reactor_set_conservation(pl_reactor *r, double relative, pl_diagnostic *diagnostic) {
    if (r == NULL || !(relative > 0.0) || !(relative < 1.0))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the conservation limit must lie between 0 and 1");
    r->conservation = relative;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_reactor_set_molar(pl_reactor *r, int species, double molar, pl_diagnostic *diagnostic) {
    if (r == NULL || species < 0 || species >= r->n)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the reactor and a species of its network are required");
    if (!isfinite(molar) || molar < 0.0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, pl_ode_time(r->ode), "a concentration must be finite and >= 0");
    if (pl_network_species_role(r->network, species) == PL_IMPLICIT)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, pl_ode_time(r->ode),
                       "%s is implicit: its concentration is folded into the rate constants",
                       pl_network_species_name(r->network, species));
    double before = r->molar[species];
    r->molar[species] = molar;
    int size = pl_ode_size(r->ode);
    double *y = malloc((size_t)size * sizeof *y);
    if (y == NULL) {
        r->molar[species] = before;
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, pl_ode_time(r->ode), "no memory to restart the reactor");
    }
    memcpy(y, pl_ode_state(r->ode), (size_t)size * sizeof *y);  /* the ledgers go on */
    for (int i = 0; i < r->dynamic; i++) y[i] = r->molar[r->index[i]];
    pl_status status = pl_ode_restart(r->ode, pl_ode_time(r->ode), y, diagnostic);
    free(y);
    if (status != PL_OK) r->molar[species] = before;
    return status;
}

double pl_reactor_time(const pl_reactor *r) { return r != NULL ? pl_ode_time(r->ode) : NAN; }

const double *pl_reactor_molar(const pl_reactor *r) { return r != NULL ? r->molar : NULL; }

void pl_reactor_get_statistics(const pl_reactor *r, pl_ode_statistics *statistics) {
    pl_ode_get_statistics(r != NULL ? r->ode : NULL, statistics);
}
