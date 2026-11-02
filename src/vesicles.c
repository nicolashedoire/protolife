/* A population of fatty acid vesicles (protolife/vesicles.h): set-up, the continuous physics and its integration. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "diagnostic_internal.h"
#include "protolife/rodas.h"
#include "vesicles_internal.h"

#define MINUTE 60.0
#define ANHYDRIDE_GONE 1e-6  /* M: the layer's face stops working below it (numerical, far below the stocks) */
#define MOST_STEPS 1000000L

pl_vesicle_constants pl_vesicle_oleate(void) {
    return (pl_vesicle_constants){.area_per_molecule = 0.32e-18, .desorption = 3.0 / MINUTE, .cvc = 0.5e-3,
                                  .rupture = 10e-3, .kelvin = 296.15, .coat_rate = 2.5, .slow_rate = 0.1,
                                  .coat_share = 0.42, .nucleation = 1.0, .born_radius = 45e-9, .critical_shear = 15.0,
                                  .division_reduced = 0.55, .least_radius = 20e-9, .apart_radius = 10e-9};
}

double pl_vesicle_molecules_for_radius(const pl_vesicle_constants *k, double radius) {
    return 2.0 * 4.0 * PL_PI * radius * radius / k->area_per_molecule;
}

pl_vesicle_shape pl_vesicle_shape_of(const pl_vesicle_constants *k, double molecules, double osmoles,
                                     double osmolarity) {
    pl_vesicle_shape shape = {.area = molecules * k->area_per_molecule / 2.0};
    shape.radius = sqrt(shape.area / (4.0 * PL_PI));
    double sphere = 4.0 / 3.0 * PL_PI * shape.radius * shape.radius * shape.radius;
    double balanced = osmolarity > 0.0 ? osmoles / osmolarity : INFINITY;
    if (balanced <= sphere) {
        shape.volume = balanced;
    } else {  /* a sphere holding more than the medium: Laplace */
        shape.volume = sphere, shape.spherical = 1;
        double pressure = PL_GAS * k->kelvin * (osmoles / sphere - osmolarity);
        shape.tension = pressure * shape.radius / 2.0;
    }
    return shape;
}

/* ---- creation ---- */

static int valid_constants(const pl_vesicle_constants *k) {
    const double positive[] = {k->area_per_molecule, k->desorption, k->cvc, k->rupture, k->kelvin, k->born_radius,
                               k->division_reduced, k->least_radius};
    const double nonnegative[] = {k->coat_rate, k->slow_rate, k->coat_share, k->nucleation, k->critical_shear,
                                  k->apart_radius, k->interface_hydrolysis, k->membrane_hydrolysis, k->uptake,
                                  k->closing};
    for (size_t i = 0; i < sizeof positive / sizeof *positive; i++)
        if (!(positive[i] > 0.0) || !isfinite(positive[i])) return 0;
    for (size_t i = 0; i < sizeof nonnegative / sizeof *nonnegative; i++)
        if (!(nonnegative[i] >= 0.0) || !isfinite(nonnegative[i])) return 0;
    return 1;
}

pl_status pl_vesicles_create(const pl_vesicle_constants *constants, int max_records, pl_vesicles **out,
                             pl_diagnostic *diagnostic) {
    if (out == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the population's address is required");
    *out = NULL;
    if (constants == NULL || !valid_constants(constants))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the vesicle constants must be finite, and positive where they size something");
    if (max_records < 0) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "max_records must be >= 0");
    pl_vesicles *p = calloc(1, sizeof *p);
    if (p == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "the population could not be allocated");
    p->k = *constants, p->max_records = max_records > 0 ? max_records : 24, p->nursery = -1, p->tolerance = 1e-8;
    p->conservation = 1e-9;
    p->record = calloc((size_t)p->max_records, sizeof *p->record);
    p->shape = calloc((size_t)p->max_records, sizeof *p->shape);
    if (p->record == NULL || p->shape == NULL) {
        pl_vesicles_free(p);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "the population could not be allocated");
    }
    *out = p;
    pl_succeed(diagnostic);
    return PL_OK;
}

void pl_vesicles_free(pl_vesicles *p) {
    if (p == NULL) return;
    for (int s = 0; s < p->solutes; s++) free(p->solute[s].name);
    free(p->solute), free(p->record), free(p->inside), free(p->shape);
    free(p->outside), free(p->leaked), free(p->feed), free(p->renewed);
    free(p->species_of), free(p->solute_of), free(p->molar), free(p->change);
    free(p->c_outside), free(p->c_inside), free(p->rate);
    pl_vesicles_free_balances(p);
    free(p);
}

/* Every per-solute array grown to `capacity`, new entries at zero (species_of at -1). */
static int grow_solutes(pl_vesicles *p, int capacity) {
    if (capacity <= p->solute_capacity) return 0;
    double **arrays[] = {&p->outside, &p->leaked, &p->feed, &p->renewed, &p->c_outside, &p->c_inside, &p->rate};
    for (size_t i = 0; i < sizeof arrays / sizeof *arrays; i++) {
        double *bigger = realloc(*arrays[i], (size_t)capacity * sizeof(double));
        if (bigger == NULL) return -1;
        for (int s = p->solute_capacity; s < capacity; s++) bigger[s] = 0.0;
        *arrays[i] = bigger;
    }
    solute_t *solutes = realloc(p->solute, (size_t)capacity * sizeof *solutes);
    if (solutes == NULL) return -1;
    p->solute = solutes;
    int *species = realloc(p->species_of, (size_t)capacity * sizeof *species);
    if (species == NULL) return -1;
    for (int s = p->solute_capacity; s < capacity; s++) species[s] = -1;
    p->species_of = species;
    p->solute_capacity = capacity;
    return 0;
}

pl_status pl_vesicles_add_solute(pl_vesicles *p, const char *name, double permeability, double osmotic, int *index,
                                 pl_diagnostic *diagnostic) {
    if (p == NULL || name == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the population and a name are required");
    if (p->frozen) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time, "solutes are added before the first vesicle");
    if (!(permeability >= 0.0) || !isfinite(permeability) || !(osmotic > 0.0) || !isfinite(osmotic))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time,
                       "the solute %s needs a permeability >= 0 and an osmotic coefficient > 0", name);
    size_t length = strlen(name) + 1;
    char *copy = malloc(length);
    int capacity = p->solutes < p->solute_capacity ? p->solute_capacity : (p->solute_capacity > 0 ? 2 * p->solute_capacity : 8);
    if (copy == NULL || grow_solutes(p, capacity) != 0) {
        free(copy);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, p->time, "no memory for the solute %s", name);
    }
    memcpy(copy, name, length);
    p->solute[p->solutes] = (solute_t){copy, permeability, osmotic};
    if (index != NULL) *index = p->solutes;
    p->solutes++;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_vesicles_attach_network(pl_vesicles *p, const pl_network *network,
                                     double (*permeability)(const pl_network *network, int species, void *context),
                                     void *context, pl_diagnostic *diagnostic) {
    if (p == NULL || network == NULL || permeability == NULL)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the population, the network and the permeability are required");
    if (p->network != NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time, "a network is already attached");
    if (p->frozen) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time, "the network is attached before the first vesicle");
    int n = pl_network_species_count(network), first = p->solutes;
    p->solute_of = malloc((size_t)(n > 0 ? n : 1) * sizeof *p->solute_of);
    p->molar = calloc((size_t)(n > 0 ? n : 1), sizeof *p->molar);
    p->change = calloc((size_t)(n > 0 ? n : 1), sizeof *p->change);
    if (p->solute_of == NULL || p->molar == NULL || p->change == NULL) {
        free(p->solute_of), free(p->molar), free(p->change);
        p->solute_of = NULL, p->molar = p->change = NULL;
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, p->time, "no memory for the network");
    }
    for (int k = 0; k < n; k++) {
        p->solute_of[k] = -1;
        if (pl_network_species_role(network, k) != PL_DYNAMIC) {
            p->molar[k] = pl_network_species_molar(network, k);
            continue;
        }
        int s;
        pl_status status = pl_vesicles_add_solute(p, pl_network_species_name(network, k), permeability(network, k, context),
                                                  1.0, &s, diagnostic);
        if (status != PL_OK) {  /* the solutes added are taken back */
            while (p->solutes > first) free(p->solute[--p->solutes].name);
            free(p->solute_of), free(p->molar), free(p->change);
            p->solute_of = NULL, p->molar = p->change = NULL;
            return status;
        }
        p->solute_of[k] = s, p->species_of[s] = k;
    }
    p->network = network;
    pl_succeed(diagnostic);
    return PL_OK;
}

int pl_vesicles_solute_of(const pl_vesicles *p, int species) {
    return p != NULL && p->network != NULL && species >= 0 && species < pl_network_species_count(p->network)
               ? p->solute_of[species] : -1;
}

int pl_vesicles_species_of(const pl_vesicles *p, int solute) {
    return p != NULL && solute >= 0 && solute < p->solutes ? p->species_of[solute] : -1;
}

/* A record appended without checks of the caller's values: its index, or -1 without room or memory. */
int pl_vesicles_append(pl_vesicles *p, double count, double molecules, const double *inside, int label) {
    if (p->records == p->max_records) return -1;
    if (p->inside == NULL) {
        p->inside = calloc((size_t)p->max_records * (size_t)(p->solutes > 0 ? p->solutes : 1), sizeof(double));
        if (p->inside == NULL) return -1;
        p->frozen = 1;
    }
    int r = p->records++;
    p->record[r] = (record_t){.count = count, .molecules = molecules, .label = label};
    double *own = pl_inside(p, r);
    for (int s = 0; s < p->solutes; s++) own[s] = inside != NULL ? inside[s] : 0.0;
    return r;
}

pl_status pl_vesicles_add(pl_vesicles *p, double count, double molecules, const double *inside, int label, int *record,
                          pl_diagnostic *diagnostic) {
    if (p == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the population is required");
    if (!(count > 0.0) || !isfinite(count) || !(molecules > 0.0) || !isfinite(molecules))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time, "a record needs a count and a number of fatty acids > 0");
    for (int s = 0; inside != NULL && s < p->solutes; s++)
        if (!(inside[s] >= 0.0) || !isfinite(inside[s]))
            return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time, "the amount of %s inside must be finite and >= 0",
                           p->solute[s].name);
    int r = pl_vesicles_append(p, count, molecules, inside, label);
    if (r < 0)
        return pl_fail(diagnostic, p->records == p->max_records ? PL_ERROR_ARGUMENT : PL_ERROR_MEMORY, 0, p->time,
                       "no room for another record (%d at most)", p->max_records);
    if (record != NULL) *record = r;
    pl_succeed(diagnostic);
    return PL_OK;
}

/* ---- the state ---- */

static int ledger_at(const pl_vesicles *p) { return p->records * (1 + p->solutes) + 1 + p->solutes; }
static int coat_at(const pl_vesicles *p) { return ledger_at(p) + p->records; }
static int intermediate_at(const pl_vesicles *p) { return coat_at(p) + p->records; }
static int pools_at(const pl_vesicles *p) { return intermediate_at(p) + p->records; }  /* forming, newborn */
static int renewed_at(const pl_vesicles *p) { return pools_at(p) + 2; }
static int anhydride_at(const pl_vesicles *p) { return renewed_at(p) + 2 * (1 + p->solutes); }
int pl_vesicles_state_size(const pl_vesicles *p) { return anhydride_at(p) + 1; }

void pl_vesicles_pack(const pl_vesicles *p, double *y) {
    int at = 0;
    for (int r = 0; r < p->records; r++) {
        const record_t *v = &p->record[r];
        const double *own = pl_inside(p, r);
        y[at++] = v->count * v->molecules / PL_AVOGADRO;
        for (int s = 0; s < p->solutes; s++) y[at++] = v->count * own[s];
    }
    y[at++] = p->monomer;
    for (int s = 0; s < p->solutes; s++) y[at++] = p->outside[s];
    for (int r = 0; r < p->records; r++) y[at++] = 0.0;
    for (int r = 0; r < p->records; r++) y[at++] = p->record[r].coat;
    for (int r = 0; r < p->records; r++) y[at++] = p->record[r].intermediate;
    y[at++] = p->forming;
    y[at++] = p->newborn;
    for (int s = 0; s < 2 * (1 + p->solutes); s++) y[at++] = 0.0;
    y[at++] = p->anhydride;
}

static void unpack(pl_vesicles *p, const double *y) {
    int at = 0;
    for (int r = 0; r < p->records; r++) {
        record_t *v = &p->record[r];
        double *own = pl_inside(p, r);
        v->molecules = y[at++] * PL_AVOGADRO / v->count;
        for (int s = 0; s < p->solutes; s++) own[s] = y[at++] / v->count;
    }
    p->monomer = y[at++];
    for (int s = 0; s < p->solutes; s++) p->outside[s] = y[at++];
    for (int r = 0; r < p->records; r++) p->record[r].released += y[at++] * PL_AVOGADRO / p->record[r].count;
    for (int r = 0; r < p->records; r++) p->record[r].coat = y[at++];
    for (int r = 0; r < p->records; r++) p->record[r].intermediate = y[at++];
    p->forming = y[at++];
    p->newborn = y[at++];
    const double *in = y + at, *out = in + 1 + p->solutes;
    p->renewed_monomer += in[0] - out[0];
    for (int s = 0; s < p->solutes; s++) p->renewed[s] += in[1 + s] - out[1 + s];
    p->anhydride = y[anhydride_at(p)];
}

/* The medium's volume depends on the vesicles' volumes, which depend on its osmolarity: a fixed point, reached in a
 * few passes since the vesicles hold a small share of the vessel. */
double pl_vesicles_read(pl_vesicles *p, const double *y, double *osmolarity_out) {
    double medium = 1.0, osmolarity = 0.0;
    const double *outside = y + p->records * (1 + p->solutes) + 1;
    for (int pass = 0; pass < 8; pass++) {
        osmolarity = 0.0;
        for (int s = 0; s < p->solutes; s++) {
            p->c_outside[s] = 1000.0 * outside[s] / medium;
            osmolarity += p->solute[s].osmotic * p->c_outside[s];
        }
        double held = 0.0;
        for (int r = 0; r < p->records; r++) {
            const double *own = y + r * (1 + p->solutes);
            double count = p->record[r].count, osmoles = 0.0;
            for (int s = 0; s < p->solutes; s++) osmoles += p->solute[s].osmotic * own[1 + s] / count;
            p->shape[r] = pl_vesicle_shape_of(&p->k, own[0] * PL_AVOGADRO / count, osmoles, osmolarity);
            held += count * p->shape[r].volume * 1000.0;
        }
        double next = 1.0 - held;
        if (fabs(next - medium) <= 1e-15) break;
        medium = next;
    }
    if (osmolarity_out != NULL) *osmolarity_out = osmolarity;
    return medium;
}

/* d[S]/dt (mol/m^3/s) of each solute from the network's reactions, at concentrations c (mol/m^3). */
static void reactions(pl_vesicles *p, const double *c, double *rate) {
    for (int s = 0; s < p->solutes; s++)
        if (p->species_of[s] >= 0) p->molar[p->species_of[s]] = c[s] / 1000.0;
    pl_network_derivative(p->network, p->molar, p->change);
    for (int s = 0; s < p->solutes; s++) rate[s] = p->species_of[s] >= 0 ? 1000.0 * p->change[p->species_of[s]] : 0.0;
}

static int derivative(double t, const double *y, double *dy, void *context) {
    (void)t;
    pl_vesicles *p = context;
    const pl_vesicle_constants *k = &p->k;
    double medium = pl_vesicles_read(p, y, NULL);
    int width = 1 + p->solutes, shared = p->records * width;
    double free_monomer = y[shared] / medium;  /* mol/L of medium */
    const double *coat = y + coat_at(p), *intermediate = y + intermediate_at(p), forming = y[pools_at(p)];
    dy[shared] = 0.0;
    for (int s = 0; s < p->solutes; s++) dy[shared + 1 + s] = 0.0;
    dy[pools_at(p)] = -k->slow_rate * fmax(forming, 0.0);
    dy[pools_at(p) + 1] = k->slow_rate * fmax(forming, 0.0);
    for (int r = 0; r < p->records; r++) {
        const double *own = y + r * width;
        double *change = dy + r * width, count = p->record[r].count;
        const pl_vesicle_shape *shape = &p->shape[r];
        /* a membrane past its rupture tension is open: it pulls no harder than at the rupture */
        double pull = fmin(shape->tension, k->rupture);
        double held = exp(-pull * k->area_per_molecule / (2.0 * PL_BOLTZMANN * k->kelvin));
        change[0] = k->desorption * own[0] * (free_monomer / k->cvc - held);
        dy[shared] -= change[0];
        dy[ledger_at(p) + r] = k->desorption * own[0] * held;
        dy[coat_at(p) + r] = -k->coat_rate * coat[r];
        dy[intermediate_at(p) + r] = -k->slow_rate * intermediate[r];
        change[0] += k->coat_rate * coat[r] + k->slow_rate * intermediate[r];
        for (int s = 0; s < p->solutes; s++) {
            double inside = shape->volume > 0.0 ? own[1 + s] / count / shape->volume : 0.0;  /* mol/m^3 */
            double flux = p->solute[s].permeability * shape->area * (p->c_outside[s] - inside);  /* mol/s each */
            change[1 + s] = count * flux;
            dy[shared + 1 + s] -= count * flux;
        }
        if (p->network != NULL && !p->record[r].inert && !p->lifeless) {  /* the reactions inside */
            for (int s = 0; s < p->solutes; s++)
                p->c_inside[s] = shape->volume > 0.0 ? own[1 + s] / count / shape->volume : 0.0;
            reactions(p, p->c_inside, p->rate);
            for (int s = 0; s < p->solutes; s++) change[1 + s] += count * shape->volume * p->rate[s];
        }
    }
    if (p->network != NULL) {  /* and in the medium */
        reactions(p, p->c_outside, p->rate);
        for (int s = 0; s < p->solutes; s++) dy[shared + 1 + s] += medium / 1000.0 * p->rate[s];
    }
    /* the anhydride layer: hydrolysed at its face (oleate to the monomer) and in the membranes (oleate staying there) */
    double anhydride = fmax(y[anhydride_at(p)], 0.0), membranes = 0.0;
    for (int r = 0; r < p->records; r++) membranes += y[r * width];
    double held = k->uptake * membranes, at_face = k->interface_hydrolysis * anhydride / (anhydride + ANHYDRIDE_GONE);
    double in_membranes = held > 0.0 ? k->membrane_hydrolysis * held * anhydride / (held + anhydride) : 0.0;
    dy[anhydride_at(p)] = p->anhydride_feed - at_face - in_membranes;
    dy[shared] += 2.0 * at_face;
    for (int r = 0; r < p->records && membranes > 0.0; r++) dy[r * width] += 2.0 * in_membranes * y[r * width] / membranes;
    /* monomer above the cvc closes into new vesicles, smoothly from the threshold */
    double excess = fmax(y[shared] - k->cvc * medium, 0.0), scale = 0.1 * k->cvc * medium;
    double above = k->closing * excess * excess / (excess + scale);
    dy[shared] -= above;
    dy[pools_at(p) + 1] += above;
    double *in = dy + renewed_at(p), *out = in + 1 + p->solutes;  /* the fresh medium in, the old out, vesicles held */
    for (int s = 0; s <= p->solutes; s++) {
        int spared = s == 0 && p->fatty_acid_held;
        in[s] = spared ? 0.0 : p->renewal * (s == 0 ? p->monomer_feed : p->feed[s - 1]) * medium;
        out[s] = spared ? 0.0 : p->renewal * y[shared + s];
        dy[shared + s] += in[s] - out[s];
    }
    return 0;
}

/* ---- integration ---- */

/* The Jacobian by differences, as Genesis's rosenbrock_differences: each component moves by 1e-6 of itself, but never
 * less than 1e-6 of a ten-thousandth of the largest component, so that one at zero beside large rates moves by more
 * than their rounding; central differences, or forward ones where going down would cross zero. */
static int jacobian(double t, const double *y, double *jac, void *context) {
    pl_vesicles *p = context;
    int n = pl_vesicles_state_size(p);
    double *shifted = p->work, *up = shifted + n, *down = shifted + 2 * n, largest = 0.0;
    memcpy(shifted, y, (size_t)n * sizeof *shifted);
    for (int j = 0; j < n; j++) largest = fmax(largest, fabs(y[j]));
    for (int j = 0; j < n; j++) {
        double h = 1e-6 * fmax(fabs(y[j]), fmax(1e-4 * largest, 1e-12));
        int central = y[j] - h >= 0.0 || y[j] < 0.0;
        shifted[j] = y[j] + h;
        derivative(t, shifted, up, p);
        shifted[j] = central ? y[j] - h : y[j];
        derivative(t, shifted, down, p);
        shifted[j] = y[j];
        for (int i = 0; i < n; i++) jac[i * n + j] = (up[i] - down[i]) / (central ? 2.0 * h : h);
    }
    return 0;
}

/* After a step: every balance changed by no more than `conservation` times the size of its terms, or the step is undone
 * and the call fails naming the balance; the changes allowed add up in the drifts. */
static pl_status check_balances(pl_vesicles *p, double dt, pl_diagnostic *diagnostic) {
    pl_vesicles_measure(p, p->after, p->scale_after);
    for (int b = 0; b < p->balances; b++) {
        double change = p->after[b] - p->before[b], size = fmax(p->scale[b], p->scale_after[b]);
        if (fabs(change) <= p->conservation * size) continue;
        const char *name = p->balance[b].name;
        double t = p->time;
        pl_vesicles_restore(p);
        return pl_fail(diagnostic, PL_ERROR_CONSERVATION, b, t,
                       "the balance of %s changed by %.3g mol/L (%.3g of its size, %.3g allowed) over the step of %g s "
                       "ending at t = %.17g s: the step is undone",
                       name, change, size > 0.0 ? fabs(change) / size : INFINITY, p->conservation, dt, t);
    }
    for (int b = 0; b < p->balances; b++)
        p->balance[b].drift += p->after[b] - p->before[b], p->balance[b].value = p->after[b];
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_vesicles_advance(pl_vesicles *p, double dt, pl_diagnostic *diagnostic) {
    if (p == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the population is required");
    if (!(dt > 0.0) || !isfinite(dt)) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time, "dt must be > 0");
    if (pl_vesicles_prepare_balances(p) != 0)
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, p->time, "no memory for the balances");
    int n = pl_vesicles_state_size(p);
    double *y = malloc(5 * (size_t)n * sizeof *y), *floor = y != NULL ? y + n : NULL;
    if (y == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, p->time, "no memory to advance the vesicles");
    p->work = y + 2 * n;
    pl_vesicles_measure(p, p->before, p->scale);
    pl_vesicles_save(p);
    pl_vesicles_pack(p, y);
    for (int i = 0; i < n; i++) floor[i] = fmax(fabs(y[i]) * 1e-6, 1e-18);
    pl_rodas_options o = {.size = n, .rhs = derivative, .jacobian = jacobian, .context = p, .tolerance = p->tolerance,
                          .floor = floor, .smallest = dt * 1e-30, .max_trials = MOST_STEPS};
    pl_status status = pl_rodas_advance(&o, p->time, y, dt, &p->step, &p->statistics, diagnostic);
    if (status == PL_OK) {
        unpack(p, y);
        p->time += dt;
        p->anhydride_fed += p->anhydride_feed * dt;
        pl_vesicles_rupture(p);
        pl_vesicles_come_apart(p);
        pl_vesicles_births(p);
        pl_vesicles_divisions(p);
        pl_vesicles_wash(p, dt);
        status = check_balances(p, dt, diagnostic);
    }
    p->work = NULL;
    free(y);
    return status;
}

/* ---- the medium and the settings ---- */

static pl_status solute_check(const pl_vesicles *p, int solute, double value, pl_diagnostic *diagnostic) {
    if (p == NULL || solute < 0 || solute >= p->solutes)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the population and one of its solutes are required");
    if (!(value >= 0.0) || !isfinite(value)) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time, "a value must be finite and >= 0");
    return PL_OK;
}

static pl_status value_check(const pl_vesicles *p, double value, pl_diagnostic *diagnostic) {
    if (p == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the population is required");
    if (!(value >= 0.0) || !isfinite(value)) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, p->time, "a value must be finite and >= 0");
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_vesicles_set_outside(pl_vesicles *p, int solute, double amount, pl_diagnostic *diagnostic) {
    pl_status status = solute_check(p, solute, amount, diagnostic);
    if (status == PL_OK) p->outside[solute] = amount, pl_succeed(diagnostic);
    return status;
}

pl_status pl_vesicles_set_feed(pl_vesicles *p, int solute, double molar, pl_diagnostic *diagnostic) {
    pl_status status = solute_check(p, solute, molar, diagnostic);
    if (status == PL_OK) p->feed[solute] = molar, pl_succeed(diagnostic);
    return status;
}

pl_status pl_vesicles_set_monomer(pl_vesicles *p, double amount, pl_diagnostic *diagnostic) {
    pl_status status = value_check(p, amount, diagnostic);
    if (status == PL_OK) p->monomer = amount;
    return status;
}

pl_status pl_vesicles_set_monomer_feed(pl_vesicles *p, double molar, pl_diagnostic *diagnostic) {
    pl_status status = value_check(p, molar, diagnostic);
    if (status == PL_OK) p->monomer_feed = molar;
    return status;
}

pl_status pl_vesicles_set_renewal(pl_vesicles *p, double rate, int fatty_acid_held, pl_diagnostic *diagnostic) {
    pl_status status = value_check(p, rate, diagnostic);
    if (status == PL_OK) p->renewal = rate, p->fatty_acid_held = fatty_acid_held != 0;
    return status;
}

pl_status pl_vesicles_set_shear(pl_vesicles *p, double rate, pl_diagnostic *diagnostic) {
    pl_status status = value_check(p, rate, diagnostic);
    if (status == PL_OK) p->shear = rate;
    return status;
}

pl_status pl_vesicles_set_washout(pl_vesicles *p, double rate, pl_diagnostic *diagnostic) {
    pl_status status = value_check(p, rate, diagnostic);
    if (status == PL_OK) p->washout = rate;
    return status;
}

pl_status pl_vesicles_set_anhydride(pl_vesicles *p, double amount, double feed, pl_diagnostic *diagnostic) {
    pl_status status = value_check(p, amount, diagnostic);
    if (status == PL_OK) status = value_check(p, feed, diagnostic);
    if (status == PL_OK) p->anhydride = amount, p->anhydride_feed = feed;
    return status;
}

pl_status pl_vesicles_set_bag(pl_vesicles *p, double litres, uint64_t seed, pl_diagnostic *diagnostic) {
    pl_status status = value_check(p, litres, diagnostic);
    if (status == PL_OK) p->bag = litres, p->draw = seed;
    return status;
}

pl_status pl_vesicles_set_tolerance(pl_vesicles *p, double relative, pl_diagnostic *diagnostic) {
    if (p == NULL || !(relative > 0.0) || !(relative < 1.0))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the tolerance must lie between 0 and 1");
    p->tolerance = relative;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_vesicles_set_conservation(pl_vesicles *p, double relative, pl_diagnostic *diagnostic) {
    if (p == NULL || !(relative > 0.0) || !(relative < 1.0))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the conservation limit must lie between 0 and 1");
    p->conservation = relative;
    pl_succeed(diagnostic);
    return PL_OK;
}

int pl_vesicles_balances(const pl_vesicles *p, const pl_balance **balances) {
    if (balances != NULL) *balances = p != NULL ? p->balance : NULL;
    return p != NULL ? p->balances : 0;
}

pl_status pl_vesicles_set_inert(pl_vesicles *p, int record, int inert, pl_diagnostic *diagnostic) {
    if (p == NULL || record < 0 || record >= p->records)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the population and one of its records are required");
    p->record[record].inert = inert != 0;
    pl_succeed(diagnostic);
    return PL_OK;
}

void pl_vesicles_set_lifeless(pl_vesicles *p, int lifeless) {
    if (p != NULL) p->lifeless = lifeless != 0;
}

/* ---- reading ---- */

double pl_vesicles_time(const pl_vesicles *p) { return p != NULL ? p->time : NAN; }
int pl_vesicles_records(const pl_vesicles *p) { return p != NULL ? p->records : 0; }
int pl_vesicles_solutes(const pl_vesicles *p) { return p != NULL ? p->solutes : 0; }
double pl_vesicles_monomer(const pl_vesicles *p) { return p != NULL ? p->monomer : NAN; }
double pl_vesicles_shear(const pl_vesicles *p) { return p != NULL ? p->shear : NAN; }

const char *pl_vesicles_solute_name(const pl_vesicles *p, int solute) {
    return p != NULL && solute >= 0 && solute < p->solutes ? p->solute[solute].name : NULL;
}

double pl_vesicles_outside(const pl_vesicles *p, int solute) {
    return p != NULL && solute >= 0 && solute < p->solutes ? p->outside[solute] : NAN;
}

pl_vesicle pl_vesicles_record(const pl_vesicles *p, int r) {
    if (p == NULL || r < 0 || r >= p->records) return (pl_vesicle){0};
    const record_t *v = &p->record[r];
    return (pl_vesicle){v->count,   v->molecules, v->released,     v->generation, v->origin,
                        v->inert,   v->label,     v->coat,         v->intermediate, pl_inside(p, r)};
}

int pl_vesicles_find(const pl_vesicles *p, int label) {
    for (int r = 0; p != NULL && r < p->records; r++)
        if (p->record[r].label == label) return r;
    return -1;
}

/* Reading needs working space; the population's own is used, since a reading changes nothing a caller sees. */
static double read_now(const pl_vesicles *p, double *osmolarity) {
    pl_vesicles *q = (pl_vesicles *)p;
    int n = pl_vesicles_state_size(p);
    double *y = malloc((size_t)n * sizeof *y);
    if (y == NULL) return NAN;
    pl_vesicles_pack(p, y);
    double medium = pl_vesicles_read(q, y, osmolarity);
    free(y);
    return medium;
}

pl_vesicle_shape pl_vesicles_shape(const pl_vesicles *p, int r) {
    if (p == NULL || r < 0 || r >= p->records || !isfinite(read_now(p, NULL))) return (pl_vesicle_shape){0};
    return p->shape[r];
}

double pl_vesicles_medium(const pl_vesicles *p, double *osmolarity) {
    if (p == NULL) return NAN;
    return read_now(p, osmolarity);
}

double pl_vesicles_membranes(const pl_vesicles *p) {
    double total = 0.0;
    for (int r = 0; p != NULL && r < p->records; r++) total += p->record[r].count * p->record[r].molecules / PL_AVOGADRO;
    return total;
}

double pl_vesicles_entering(const pl_vesicles *p) {
    if (p == NULL) return NAN;
    double total = p->forming + p->newborn;
    for (int r = 0; r < p->records; r++) total += p->record[r].coat + p->record[r].intermediate;
    return total;
}

double pl_vesicles_forming(const pl_vesicles *p) { return p != NULL ? p->forming + p->newborn : NAN; }

double pl_vesicles_fatty_acid(const pl_vesicles *p) {
    return p != NULL ? p->monomer + pl_vesicles_entering(p) + pl_vesicles_membranes(p) : NAN;
}

double pl_vesicles_solute_total(const pl_vesicles *p, int solute) {
    if (p == NULL || solute < 0 || solute >= p->solutes) return NAN;
    double total = p->outside[solute];
    for (int r = 0; r < p->records; r++) total += p->record[r].count * pl_inside(p, r)[solute];
    return total;
}

pl_vesicle_ledger pl_vesicles_ledger(const pl_vesicles *p) {
    if (p == NULL) return (pl_vesicle_ledger){0};
    return (pl_vesicle_ledger){p->ruptures, p->feeds,           p->divided, p->born,   p->apart,
                               p->washed,   p->renewed_monomer, p->fed,     p->renewed, p->leaked};
}

void pl_vesicles_get_statistics(const pl_vesicles *p, pl_ode_statistics *statistics) {
    if (statistics == NULL) return;
    *statistics = p != NULL ? p->statistics : (pl_ode_statistics){0};
}
