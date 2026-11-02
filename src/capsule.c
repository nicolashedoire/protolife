/* One protocell after Taneja and Higgs 2025 (protolife/capsule.h). */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "capsule_internal.h"
#include "diagnostic_internal.h"

#define PI 3.14159265358979323846
#define AVOGADRO 6.02e23  /* the value of the authors' program, for their Q = 0.0249 */

struct pl_capsule {
    pl_capsule_parameters k;
    const pl_network *network;
    int dynamic, *index;        /* the network's dynamic species, inside the cell */
    double *mu, *outside;       /* per dynamic species */
    double *molar, *change;     /* per network species: working space of the reactions */
    double *y, *stage, *slope[4];  /* V, A+, A-, each inside concentration, the lipid in solution */
    double DA0, Q, time;
    int keep;                   /* the daughter followed next: 1 the larger, -1 the smaller */
    uint64_t random;
    pl_capsule_event *event;
    int events, event_capacity;
};

pl_capsule_parameters pl_capsule_taneja_higgs(void) {
    return (pl_capsule_parameters){.lambda = 1.0, .rmem = 0.1, .rflip = 0.1, .rdiv = 10.0, .cstar = 0.01,
                                   .elip = 1.01 * 0.01, .R0 = 1e-6, .d = 1e-9, .a = 2.0e-19, .kappa = 1.0,
                                   .beta = 0.375, .burst_area = 1.1, .burst_kept = 0.9, .largest = 2000.0};
}

/* splitmix64, then 53 bits: a uniform number in [0, 1), the same on every platform. */
static double uniform(pl_capsule *c) {
    uint64_t z = (c->random += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return (double)((z ^ (z >> 31)) >> 11) * 0x1.0p-53;
}

void pl_capsule_free(pl_capsule *c) {
    if (c == NULL) return;
    free(c->index), free(c->mu), free(c->outside), free(c->molar), free(c->change), free(c->y), free(c->event);
    free(c);
}

static int valid(const pl_capsule_parameters *k) {
    const double positive[] = {k->cstar, k->R0, k->d, k->a, k->kappa, k->burst_area, k->burst_kept, k->largest};
    const double nonnegative[] = {k->lambda, k->rmem, k->rflip, k->rdiv, k->elip, k->beta, k->kappa_gauss};
    for (size_t i = 0; i < sizeof positive / sizeof *positive; i++)
        if (!(positive[i] > 0.0) || !isfinite(positive[i])) return 0;
    for (size_t i = 0; i < sizeof nonnegative / sizeof *nonnegative; i++)
        if (!(nonnegative[i] >= 0.0) || !isfinite(nonnegative[i])) return 0;
    return k->burst_area > 1.0 && k->burst_kept < 1.0;
}

pl_status pl_capsule_create(const pl_capsule_parameters *k, const pl_network *network, const double *permeation,
                            const double *outside, uint64_t seed, pl_capsule **out, pl_diagnostic *diagnostic) {
    if (out == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the capsule's address is required");
    *out = NULL;
    if (k == NULL || !valid(k) || network == NULL || permeation == NULL || outside == NULL)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN,
                       "a capsule needs valid parameters, a network, and its permeations and outside concentrations");
    int n = pl_network_species_count(network);
    pl_capsule *c = calloc(1, sizeof *c);
    if (c == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "the capsule could not be allocated");
    c->k = *k, c->network = network, c->keep = -1, c->random = seed;
    c->DA0 = 4.0 * k->d / k->R0;
    c->Q = (4.0 * PI * k->R0 * k->R0 / k->a) / (AVOGADRO * 4.0 / 3.0 * PI * pow(k->R0, 3.0) * 1000.0);
    c->index = malloc((size_t)(n > 0 ? n : 1) * sizeof *c->index);
    c->molar = malloc((size_t)(n > 0 ? n : 1) * sizeof *c->molar), c->change = malloc((size_t)(n > 0 ? n : 1) * sizeof *c->change);
    if (c->index == NULL || c->molar == NULL || c->change == NULL) return pl_capsule_free(c), pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory");
    for (int s = 0; s < n; s++) {
        c->molar[s] = pl_network_species_molar(network, s);
        if (pl_network_species_role(network, s) == PL_DYNAMIC) c->index[c->dynamic++] = s;
    }
    int size = 4 + c->dynamic;
    c->mu = malloc((size_t)(c->dynamic > 0 ? c->dynamic : 1) * sizeof *c->mu);
    c->outside = malloc((size_t)(c->dynamic > 0 ? c->dynamic : 1) * sizeof *c->outside);
    c->y = malloc(6 * (size_t)size * sizeof *c->y);
    if (c->mu == NULL || c->outside == NULL || c->y == NULL) return pl_capsule_free(c), pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory");
    c->stage = c->y + size;
    for (int j = 0; j < 4; j++) c->slope[j] = c->y + (size_t)(2 + j) * (size_t)size;
    for (int i = 0; i < c->dynamic; i++) {
        c->mu[i] = permeation[c->index[i]], c->outside[i] = outside[c->index[i]];
        if (!(c->mu[i] >= 0.0) || !(c->outside[i] >= 0.0))
            return pl_capsule_free(c), pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "permeations and outside concentrations must be >= 0");
        c->y[3 + i] = c->molar[c->index[i]];
    }
    /* a sphere of radius R0 at rest: leaflets as its curvature asks, C* of lipid inside */
    c->y[0] = 1.0, c->y[1] = 1.0 + 0.5 * c->DA0, c->y[2] = 1.0 - 0.5 * c->DA0, c->y[3 + c->dynamic] = k->cstar;
    *out = c;
    pl_succeed(diagnostic);
    return PL_OK;
}

/* dy/dt. With `fixed`, the chemistry alone in a compartment of A/V = A0/V0 (the authors' equilibration). */
static void derivative(pl_capsule *c, const double *y, double *dy, int fixed) {
    const pl_capsule_parameters *k = &c->k;
    int n = c->dynamic;
    for (int i = 0; i < n; i++) c->molar[c->index[i]] = y[3 + i];
    pl_network_derivative(c->network, c->molar, c->change);
    if (fixed) {
        dy[0] = dy[1] = dy[2] = dy[3 + n] = 0.0;
        for (int i = 0; i < n; i++) dy[3 + i] = c->change[c->index[i]] + c->mu[i] * (c->outside[i] - y[3 + i]);
        return;
    }
    double volume = y[0], lipid = y[3 + n];
    capsule_form f = pl_capsule_form(volume, 0.5 * (y[1] + y[2]));
    double asked = (f.radius + 0.25 * f.length) * c->DA0, outer = f.actual + 0.5 * asked, inner = f.actual - 0.5 * asked;
    double excess = lipid - k->elip;
    for (int i = 0; i < n; i++) excess += y[3 + i] - c->outside[i];
    double growth = k->lambda * f.actual * excess, flip = k->rflip * (y[1] * inner - y[2] * outer) / f.actual;
    dy[0] = growth;
    dy[1] = k->rmem * (outer * k->elip / k->cstar - y[1]) - flip;
    dy[2] = k->rmem * (inner * lipid / k->cstar - y[2]) + flip;
    for (int i = 0; i < n; i++)
        dy[3 + i] = c->change[c->index[i]] + c->mu[i] * (f.actual / volume) * (c->outside[i] - y[3 + i]) -
                    y[3 + i] / volume * growth;
    dy[3 + n] = -(c->Q * k->rmem / volume) * (inner * lipid / k->cstar - y[2]) - lipid / volume * growth;
}

/* One classical Runge-Kutta step. (The authors' program evaluates the fourth stage at half a step; with their small
 * steps both converge to the same solution.) */
static void runge_kutta(pl_capsule *c, double dt, int fixed) {
    int size = 4 + c->dynamic;
    const double weight[3] = {0.5, 0.5, 1.0};
    derivative(c, c->y, c->slope[0], fixed);
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < size; i++) c->stage[i] = c->y[i] + weight[j] * dt * c->slope[j][i];
        derivative(c, c->stage, c->slope[j + 1], fixed);
    }
    for (int i = 0; i < size; i++)
        c->y[i] += dt * (c->slope[0][i] + 2.0 * c->slope[1][i] + 2.0 * c->slope[2][i] + c->slope[3][i]) / 6.0;
}

static int record(pl_capsule *c, pl_capsule_event e) {
    if (c->events == c->event_capacity) {
        int capacity = c->event_capacity > 0 ? 2 * c->event_capacity : 64;
        pl_capsule_event *bigger = realloc(c->event, (size_t)capacity * sizeof *bigger);
        if (bigger == NULL) return -1;
        c->event = bigger, c->event_capacity = capacity;
    }
    c->event[c->events++] = e;
    return 0;
}

/* The energy rule: a sphere of a random share x < 1/2 of the volume and the rest, the area difference shared in
 * proportion to area; done when the daughters' energy is below the parent's. 1 when it divided, -1 for a thin tube. */
static int try_division(pl_capsule *c, double reduced) {
    double V = c->y[0], outer = c->y[1], inner = c->y[2], lipid = 0.5 * (outer + inner);
    capsule_form f = pl_capsule_form(V, lipid);
    if (f.radius < 0.01) return -1;  /* a very long thin tube: the authors stop here */
    double parent = pl_capsule_energy(&c->k, c->DA0, V, outer, inner), x = 0.5 * uniform(c);
    double V1 = V * x, V2 = V - V1, A1 = pow(V1, 2.0 / 3.0), A2 = lipid - A1, difference = outer - inner;
    if (!(A2 >= pow(V2, 2.0 / 3.0))) return 0;
    double n1 = difference * A1 / lipid, n2 = difference * A2 / lipid;
    double E1 = pl_capsule_energy(&c->k, c->DA0, V1, A1 + 0.5 * n1, A1 - 0.5 * n1);
    double E2 = pl_capsule_energy(&c->k, c->DA0, V2, A2 + 0.5 * n2, A2 - 0.5 * n2);
    if (!(E1 + E2 + 4.0 * PI * c->k.kappa_gauss - parent < 0.0)) return 0;
    c->keep = -c->keep;
    double share = c->keep > 0 ? V2 / V : V1 / V, A = c->keep > 0 ? A2 : A1, n = c->keep > 0 ? n2 : n1;
    record(c, (pl_capsule_event){PL_CAPSULE_DIVISION, c->time, reduced, V, share});
    c->y[0] = share * V, c->y[1] = A + 0.5 * n, c->y[2] = A - 0.5 * n;
    return 1;
}

static double reduced_volume(const pl_capsule *c) { return c->y[0] / pow(0.5 * (c->y[1] + c->y[2]), 1.5); }

pl_status pl_capsule_advance(pl_capsule *c, double duration, double dt, pl_diagnostic *diagnostic) {
    if (c == NULL || !(dt > 0.0) || !(duration >= 0.0))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the capsule, a step > 0 and a duration >= 0 are required");
    long steps = lround(duration / dt);
    double burst = pow(c->k.burst_area, 1.5);
    for (long step = 0; step < steps; step++) {
        runge_kutta(c, dt, 0);
        c->time += dt;
        double v = reduced_volume(c);
        if (c->k.equal_division && v <= 0.5) {
            record(c, (pl_capsule_event){PL_CAPSULE_DIVISION, c->time, v, c->y[0], 0.5});
            c->y[0] *= 0.5, c->y[1] *= 0.5, c->y[2] *= 0.5;
        } else if (!c->k.equal_division && v < 1.0 && uniform(c) < c->k.rdiv * dt) {
            if (try_division(c, v) < 0)
                return pl_fail(diagnostic, PL_ERROR_MODEL, 0, c->time, "the cell became a tube thinner than R0/100 at t = %g", c->time);
            v = reduced_volume(c);
        }
        if (v > burst) {
            record(c, (pl_capsule_event){PL_CAPSULE_BURST, c->time, v, c->y[0], c->k.burst_kept});
            c->y[0] *= c->k.burst_kept;
        }
        if (c->y[0] > c->k.largest)
            return pl_fail(diagnostic, PL_ERROR_MODEL, 0, c->time, "the cell grew past %g V0 at t = %g", c->k.largest, c->time);
        for (int i = 0; i < 4 + c->dynamic; i++)
            if (!isfinite(c->y[i])) return pl_fail(diagnostic, PL_ERROR_MODEL, 0, c->time, "the state is not finite at t = %g", c->time);
    }
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_capsule_settle(pl_capsule *c, double duration, double dt, pl_diagnostic *diagnostic) {
    if (c == NULL || !(dt > 0.0) || !(duration >= 0.0))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the capsule, a step > 0 and a duration >= 0 are required");
    long steps = lround(duration / dt);
    for (long step = 0; step < steps; step++) runge_kutta(c, dt, 1);
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_capsule_set_molar(pl_capsule *c, int species, double molar, pl_diagnostic *diagnostic) {
    for (int i = 0; c != NULL && i < c->dynamic; i++)
        if (c->index[i] == species) {
            if (!(molar >= 0.0) || !isfinite(molar))
                return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, c->time, "a concentration must be finite and >= 0");
            c->y[3 + i] = molar;
            pl_succeed(diagnostic);
            return PL_OK;
        }
    return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the capsule and one of its network's dynamic species are required");
}

double pl_capsule_time(const pl_capsule *c) { return c != NULL ? c->time : NAN; }

pl_capsule_shape pl_capsule_get_shape(const pl_capsule *c) {
    if (c == NULL) return (pl_capsule_shape){0};
    double V = c->y[0], lipid = 0.5 * (c->y[1] + c->y[2]);
    capsule_form f = pl_capsule_form(V, lipid);
    return (pl_capsule_shape){V, f.actual, lipid, c->y[1] - c->y[2], f.radius, f.length, V / pow(lipid, 1.5), lipid / V,
                              pl_capsule_energy(&c->k, c->DA0, V, c->y[1], c->y[2]), f.swollen};
}

double pl_capsule_molar(const pl_capsule *c, int species) {
    for (int i = 0; c != NULL && i < c->dynamic; i++)
        if (c->index[i] == species) return c->y[3 + i];
    return c != NULL && species >= 0 && species < pl_network_species_count(c->network) ? c->molar[species] : NAN;
}

double pl_capsule_lipid(const pl_capsule *c) { return c != NULL ? c->y[3 + c->dynamic] : NAN; }

int pl_capsule_events(const pl_capsule *c, const pl_capsule_event **events) {
    if (events != NULL) *events = c != NULL ? c->event : NULL;
    return c != NULL ? c->events : 0;
}
