/* A well-stirred reactor over pl_ode (protolife/reactor.h), after Genesis's chemistry/mean_field.c. */
#include "protolife/reactor.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "diagnostic_internal.h"

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
    return 0;
}

static int jacobian(double t, const double *y, double *out, void *context) {
    (void)t;
    pl_reactor *r = context;
    expand(r, y);
    pl_network_jacobian(r->network, r->scratch, r->jacobian);
    for (int i = 0; i < r->dynamic; i++)
        for (int j = 0; j < r->dynamic; j++)
            out[i * r->dynamic + j] = r->jacobian[r->index[i] * r->n + r->index[j]] - (i == j ? r->dilution : 0.0);
    return 0;
}

void pl_reactor_free(pl_reactor *r) {
    if (r == NULL) return;
    pl_ode_free(r->ode);
    free(r->index), free(r->feed), free(r->molar), free(r->scratch), free(r->derivative), free(r->jacobian);
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
    for (int i = 0; i < r->dynamic; i++)
        y0[i] = r->molar[r->index[i]], absolute[i] = o.absolute_tolerance > 0.0 ? o.absolute_tolerance : 1e-14;
    if (r->dynamic == 0) {
        free(y0), free(absolute), pl_reactor_free(r);
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, 0.0, "the network has no dynamic species to integrate");
    }
    {
        pl_ode_options ode = {.size = r->dynamic, .rhs = rhs, .jacobian = jacobian, .context = r,
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
    pl_status status = pl_ode_advance(r->ode, t_end, NULL, diagnostic);
    if (status != PL_OK) return status;
    const double *y = pl_ode_state(r->ode);
    for (int i = 0; i < r->dynamic; i++) r->molar[r->index[i]] = y[i];
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
    double *y = malloc((size_t)r->dynamic * sizeof *y);
    if (y == NULL) {
        r->molar[species] = before;
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, pl_ode_time(r->ode), "no memory to restart the reactor");
    }
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
