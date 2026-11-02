/* Reaction networks (protolife/network.h), after Genesis's chemistry/network.c and the rate law of mean_field.c. */
#include "protolife/network.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diagnostic_internal.h"
#include "formula_internal.h"

enum { LEFT, RIGHT, IMPLICIT_LEFT, IMPLICIT_RIGHT, SIDES };

typedef struct {
    int species, coefficient;
} term;

typedef struct {
    char *name, *formula;
    int charge;
    pl_role role;
    double molar;
    int atoms[PL_ELEMENTS];
} species;

typedef struct {
    term *terms;        /* the four sides one after the other */
    int count[SIDES];
    double rate;
    int lumped;
    char *source, *equation;
} reaction;

struct pl_network {
    species *species;
    int species_count, species_capacity;
    reaction *reactions;
    int reaction_count, reaction_capacity;
    double temperature;
};

static char *copy_text(const char *text) {
    size_t n = strlen(text) + 1;
    char *copy = malloc(n);
    if (copy != NULL) memcpy(copy, text, n);
    return copy;
}

pl_status pl_network_create(pl_network **network, pl_diagnostic *diagnostic) {
    if (network == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the network's address is required");
    *network = calloc(1, sizeof **network);
    if (*network == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "the network could not be allocated");
    (*network)->temperature = 298.15;
    pl_succeed(diagnostic);
    return PL_OK;
}

void pl_network_free(pl_network *network) {
    if (network == NULL) return;
    for (int s = 0; s < network->species_count; s++) free(network->species[s].name), free(network->species[s].formula);
    for (int r = 0; r < network->reaction_count; r++) {
        reaction *x = &network->reactions[r];
        free(x->terms), free(x->source), free(x->equation);
    }
    free(network->species);
    free(network->reactions);
    free(network);
}

int pl_network_find(const pl_network *network, const char *name) {
    if (network == NULL || name == NULL) return -1;
    for (int s = 0; s < network->species_count; s++)
        if (strcmp(network->species[s].name, name) == 0) return s;
    return -1;
}

static int valid_name(const char *name) {
    if (name == NULL || *name == '\0' || isdigit((unsigned char)*name) || strcmp(name, "+") == 0) return 0;
    for (const char *c = name; *c != '\0'; c++)
        if (isspace((unsigned char)*c) || *c == ';' || *c == '#') return 0;
    return strstr(name, "->") == NULL;
}

/* Makes room for one more element in an array of `size`-byte elements. */
static int grow(void **array, int *capacity, int count, size_t size) {
    if (count < *capacity) return 0;
    int next = *capacity > 0 ? 2 * *capacity : 16;
    void *bigger = realloc(*array, (size_t)next * size);
    if (bigger == NULL) return -1;
    *array = bigger, *capacity = next;
    return 0;
}

pl_status pl_network_add_species(pl_network *network, const char *name, const char *formula, int charge, pl_role role,
                                 double molar, pl_diagnostic *diagnostic) {
    if (network == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the network is required");
    if (!valid_name(name))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN,
                       "a species name must be non-empty, without blanks, ';', '#' or '->', and not start with a digit");
    if (pl_network_find(network, name) >= 0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the species %s already exists", name);
    if (role != PL_DYNAMIC && role != PL_RESERVOIR && role != PL_IMPLICIT)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the role of %s is not a pl_role", name);
    if (!isfinite(molar) || molar < 0.0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the concentration of %s must be finite and >= 0", name);
    species s = {.charge = charge, .role = role, .molar = molar};
    if (pl_formula_parse(formula, s.atoms) != 0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the formula of %s (%s) cannot be read", name,
                       formula != NULL ? formula : "none");
    s.name = copy_text(name), s.formula = copy_text(formula);
    if (s.name == NULL || s.formula == NULL ||
        grow((void **)&network->species, &network->species_capacity, network->species_count, sizeof s) != 0) {
        free(s.name), free(s.formula);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory for the species %s", name);
    }
    network->species[network->species_count++] = s;
    pl_succeed(diagnostic);
    return PL_OK;
}

/* The terms of one side, "2 X + Y" or "2X + Y" or empty, into terms[] (at most `room`): their number, or -1. */
static int parse_side(const pl_network *network, const char *text, term *terms, int room, char *problem,
                      size_t problem_size) {
    int count = 0;
    const char *p = text;
    for (;;) {
        while (isspace((unsigned char)*p)) p++;
        if (*p == '\0') break;
        if (count > 0) {
            if (*p != '+') return snprintf(problem, problem_size, "expected '+' before \"%s\"", p), -1;
            p++;
            while (isspace((unsigned char)*p)) p++;
        }
        long coefficient = 1;
        if (isdigit((unsigned char)*p)) {
            char *end;
            errno = 0;
            coefficient = strtol(p, &end, 10);
            if (errno != 0 || coefficient <= 0 || coefficient > 1000)
                return snprintf(problem, problem_size, "a coefficient must be between 1 and 1000"), -1;
            p = end;
            while (isspace((unsigned char)*p)) p++;
        }
        const char *start = p;
        while (*p != '\0' && !isspace((unsigned char)*p)) p++;
        char name[256];
        size_t length = (size_t)(p - start);
        if (length == 0 || length >= sizeof name) return snprintf(problem, problem_size, "a term has no name"), -1;
        memcpy(name, start, length), name[length] = '\0';
        int s = pl_network_find(network, name);
        if (s < 0) return snprintf(problem, problem_size, "unknown species %.200s", name), -1;
        if (count == room) return snprintf(problem, problem_size, "too many terms"), -1;
        terms[count++] = (term){s, (int)coefficient};
    }
    return count;
}

/* "left -> right" into two sides of `terms` (room for `room` in all): 0, or -1 with the problem. */
static int parse_equation(const pl_network *network, const char *equation, term *terms, int room, int count[2],
                          char *problem, size_t problem_size) {
    const char *arrow = strstr(equation, "->");
    if (arrow == NULL || strstr(arrow + 2, "->") != NULL)
        return snprintf(problem, problem_size, "an equation needs exactly one \"->\""), -1;
    size_t length = (size_t)(arrow - equation);
    char *left = malloc(length + 1);
    if (left == NULL) return snprintf(problem, problem_size, "no memory"), -1;
    memcpy(left, equation, length), left[length] = '\0';
    count[0] = parse_side(network, left, terms, room, problem, problem_size);
    free(left);
    if (count[0] < 0) return -1;
    count[1] = parse_side(network, arrow + 2, terms + count[0], room - count[0], problem, problem_size);
    return count[1] < 0 ? -1 : 0;
}

/* Terms an equation may hold: a bound from its length, since every term takes at least one character. */
static int term_room(const char *equation) { return equation != NULL ? (int)strlen(equation) / 2 + 2 : 0; }

pl_status pl_network_add_reaction(pl_network *network, const char *equation, double rate, const char *implicit,
                                  int lumped, const char *source, pl_diagnostic *diagnostic) {
    if (network == NULL || equation == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the network and the equation are required");
    if (!isfinite(rate) || rate < 0.0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the rate constant of %s must be finite and >= 0", equation);
    if (source == NULL || *source == '\0')
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the reaction %s has no source", equation);
    int room = term_room(equation) + term_room(implicit);
    reaction x = {.rate = rate, .lumped = lumped != 0, .terms = malloc((size_t)room * sizeof(term))};
    if (x.terms == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory for the reaction %s", equation);
    char problem[256] = "";
    int used = 0;
    if (parse_equation(network, equation, x.terms, room, &x.count[LEFT], problem, sizeof problem) == 0) {
        used = x.count[LEFT] + x.count[RIGHT];
        if (implicit != NULL &&
            parse_equation(network, implicit, x.terms + used, room - used, &x.count[IMPLICIT_LEFT], problem,
                           sizeof problem) != 0)
            used = -1;
    } else {
        used = -1;
    }
    if (used < 0) {
        free(x.terms);
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "cannot read the reaction %s%s%s: %s", equation,
                       implicit != NULL ? " ; implicit=" : "", implicit != NULL ? implicit : "", problem);
    }
    x.source = copy_text(source), x.equation = copy_text(equation);
    if (x.source == NULL || x.equation == NULL ||
        grow((void **)&network->reactions, &network->reaction_capacity, network->reaction_count, sizeof x) != 0) {
        free(x.terms), free(x.source), free(x.equation);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory for the reaction %s", equation);
    }
    network->reactions[network->reaction_count++] = x;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_network_set_temperature(pl_network *network, double kelvin, pl_diagnostic *diagnostic) {
    if (network == NULL || !isfinite(kelvin) || kelvin <= 0.0)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the temperature must be positive and finite");
    network->temperature = kelvin;
    pl_succeed(diagnostic);
    return PL_OK;
}

/* ---- the text format ---- */

static char *trim(char *s) {
    while (isspace((unsigned char)*s)) s++;
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n - 1])) s[--n] = '\0';
    return s;
}

static pl_status line_species(pl_network *network, char *rest, pl_diagnostic *diagnostic) {
    char name[256], formula[256], role[32], extra[2];
    int charge;
    double molar;
    if (sscanf(rest, "%255s %255s %d %31s %lf %1s", name, formula, &charge, role, &molar, extra) != 5)
        return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "expected: species <name> <formula> <charge> <role> <mol/L>");
    pl_role r = strcmp(role, "dynamic") == 0 ? PL_DYNAMIC : strcmp(role, "reservoir") == 0 ? PL_RESERVOIR
                : strcmp(role, "implicit") == 0 ? PL_IMPLICIT : (pl_role)-1;
    if (r == (pl_role)-1)
        return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "the role %s is not dynamic, reservoir or implicit", role);
    return pl_network_add_species(network, name, formula, charge, r, molar, diagnostic);
}

static pl_status line_reaction(pl_network *network, char *rest, pl_diagnostic *diagnostic) {
    char *field = strchr(rest, ';'), *equation = rest, *implicit = NULL, *source = NULL;
    double rate = NAN;
    int lumped = 0;
    if (field != NULL) *field++ = '\0';
    while (field != NULL) {
        char *next = strchr(field, ';');
        if (next != NULL) *next++ = '\0';
        char *f = trim(field), *end = NULL;
        if (strncmp(f, "k=", 2) == 0) {
            rate = strtod(f + 2, &end);
            if (end == f + 2 || *trim(end) != '\0')
                return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "the rate constant \"%s\" is not a number", f + 2);
        } else if (strncmp(f, "implicit=", 9) == 0) {
            implicit = f + 9;
        } else if (strcmp(f, "lumped") == 0) {
            lumped = 1;
        } else if (strncmp(f, "source=", 7) == 0) {
            source = f + 7;
        } else {
            return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "unknown field \"%s\"", f);
        }
        field = next;
    }
    if (isnan(rate)) return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "the reaction has no k=");
    return pl_network_add_reaction(network, trim(equation), rate, implicit, lumped, source, diagnostic);
}

static pl_status parse_line(pl_network *network, char *text, pl_diagnostic *diagnostic) {
    text = trim(text);
    if (*text == '\0' || *text == '#') return PL_OK;
    if (strncmp(text, "species", 7) == 0 && isspace((unsigned char)text[7])) return line_species(network, text + 8, diagnostic);
    if (strncmp(text, "reaction", 8) == 0 && isspace((unsigned char)text[8])) return line_reaction(network, text + 9, diagnostic);
    if (strncmp(text, "temperature", 11) == 0 && isspace((unsigned char)text[11])) {
        char *end;
        double kelvin = strtod(text + 12, &end);
        if (end == text + 12 || *trim(end) != '\0')
            return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "the temperature is not a number");
        return pl_network_set_temperature(network, kelvin, diagnostic);
    }
    return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "a line starts with species, reaction or temperature");
}

pl_status pl_network_parse(const char *text, pl_network **out, pl_diagnostic *diagnostic) {
    if (out == NULL || text == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the text and the network's address are required");
    *out = NULL;
    pl_network *network;
    pl_status status = pl_network_create(&network, diagnostic);
    if (status != PL_OK) return status;
    char *copy = copy_text(text);
    if (copy == NULL) {
        pl_network_free(network);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory to read the network");
    }
    int number = 1;
    for (char *line = copy, *next; line != NULL && status == PL_OK; line = next, number++) {
        next = strchr(line, '\n');
        if (next != NULL) *next++ = '\0';
        pl_diagnostic local;
        status = parse_line(network, line, &local);
        if (status != PL_OK) {
            status = local.status == PL_ERROR_MEMORY ? PL_ERROR_MEMORY : PL_ERROR_INPUT;
            pl_fail(diagnostic, status, number, NAN, "line %d: %s", number, local.message);
        }
    }
    free(copy);
    if (status != PL_OK) {
        pl_network_free(network);
        return status;
    }
    *out = network;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_network_load(const char *path, pl_network **out, pl_diagnostic *diagnostic) {
    if (out == NULL || path == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the path and the network's address are required");
    *out = NULL;
    FILE *in = fopen(path, "rb");
    if (in == NULL) return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "cannot open %s", path);
    char *text = NULL;
    size_t length = 0, capacity = 0;
    for (;;) {
        if (length + 4096 + 1 > capacity) {
            capacity = 2 * capacity + 8192;
            char *bigger = realloc(text, capacity);
            if (bigger == NULL) {
                free(text), fclose(in);
                return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory to read %s", path);
            }
            text = bigger;
        }
        size_t got = fread(text + length, 1, 4096, in);
        length += got;
        if (got < 4096) break;
    }
    int failed = ferror(in);
    fclose(in);
    if (failed) {
        free(text);
        return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "cannot read %s", path);
    }
    text[length] = '\0';
    pl_status status = pl_network_parse(text, out, diagnostic);
    free(text);
    if (status != PL_OK && diagnostic != NULL) {
        char message[PL_MESSAGE_SIZE];
        snprintf(message, sizeof message, "%s", diagnostic->message);
        pl_fail(diagnostic, status, diagnostic->code, NAN, "%s, %s", path, message);
    }
    return status;
}

/* ---- the balance ---- */

static void add_side(const pl_network *network, const reaction *x, int side, int atoms[PL_ELEMENTS], int *charge) {
    int first = 0;
    for (int k = 0; k < side; k++) first += x->count[k];
    for (int t = first; t < first + x->count[side]; t++) {
        const species *s = &network->species[x->terms[t].species];
        for (int e = 0; e < PL_ELEMENTS; e++) atoms[e] += x->terms[t].coefficient * s->atoms[e];
        *charge += x->terms[t].coefficient * s->charge;
    }
}

pl_status pl_network_check_balance(const pl_network *network, pl_diagnostic *diagnostic) {
    if (network == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "the network is required");
    for (int r = 0; r < network->reaction_count; r++) {
        const reaction *x = &network->reactions[r];
        if (x->lumped) continue;
        int left[PL_ELEMENTS] = {0}, right[PL_ELEMENTS] = {0}, left_charge = 0, right_charge = 0;
        add_side(network, x, LEFT, left, &left_charge), add_side(network, x, IMPLICIT_LEFT, left, &left_charge);
        add_side(network, x, RIGHT, right, &right_charge), add_side(network, x, IMPLICIT_RIGHT, right, &right_charge);
        for (int e = 0; e < PL_ELEMENTS; e++)
            if (left[e] != right[e])
                return pl_fail(diagnostic, PL_ERROR_CONSERVATION, r, NAN,
                               "reaction %d (%s; %s) is not balanced: %d %s on the left, %d on the right", r + 1,
                               x->equation, x->source, left[e], pl_element_symbol(e), right[e]);
        if (left_charge != right_charge)
            return pl_fail(diagnostic, PL_ERROR_CONSERVATION, r, NAN,
                           "reaction %d (%s; %s) is not balanced: charge %d on the left, %d on the right", r + 1,
                           x->equation, x->source, left_charge, right_charge);
    }
    pl_succeed(diagnostic);
    return PL_OK;
}

int pl_network_conserves(const pl_network *network, const char *element) {
    int e = pl_element_find(element);
    if (network == NULL || e < 0) return 0;
    for (int s = 0; s < network->species_count; s++)
        if (network->species[s].role != PL_DYNAMIC && network->species[s].atoms[e] != 0) return 0;
    for (int r = 0; r < network->reaction_count; r++) {
        const reaction *x = &network->reactions[r];
        long left = 0, right = 0;
        for (int t = 0; t < x->count[LEFT] + x->count[RIGHT]; t++) {
            long n = (long)x->terms[t].coefficient * network->species[x->terms[t].species].atoms[e];
            if (t < x->count[LEFT]) left += n;
            else right += n;
        }
        if (left != right) return 0;
    }
    return 1;
}

/* ---- accessors ---- */

int pl_network_species_count(const pl_network *network) { return network != NULL ? network->species_count : 0; }
int pl_network_reaction_count(const pl_network *network) { return network != NULL ? network->reaction_count : 0; }
double pl_network_temperature(const pl_network *network) { return network != NULL ? network->temperature : NAN; }

static const species *species_at(const pl_network *network, int s) {
    return network != NULL && s >= 0 && s < network->species_count ? &network->species[s] : NULL;
}

static const reaction *reaction_at(const pl_network *network, int r) {
    return network != NULL && r >= 0 && r < network->reaction_count ? &network->reactions[r] : NULL;
}

const char *pl_network_species_name(const pl_network *network, int s) {
    const species *x = species_at(network, s);
    return x != NULL ? x->name : NULL;
}

const char *pl_network_species_formula(const pl_network *network, int s) {
    const species *x = species_at(network, s);
    return x != NULL ? x->formula : NULL;
}

int pl_network_species_charge(const pl_network *network, int s) {
    const species *x = species_at(network, s);
    return x != NULL ? x->charge : 0;
}

int pl_network_species_atoms(const pl_network *network, int s, const char *element) {
    const species *x = species_at(network, s);
    int e = pl_element_find(element);
    return x != NULL && e >= 0 ? x->atoms[e] : 0;
}

pl_role pl_network_species_role(const pl_network *network, int s) {
    const species *x = species_at(network, s);
    return x != NULL ? x->role : PL_DYNAMIC;
}

double pl_network_species_molar(const pl_network *network, int s) {
    const species *x = species_at(network, s);
    return x != NULL ? x->molar : NAN;
}

double pl_network_reaction_constant(const pl_network *network, int r) {
    const reaction *x = reaction_at(network, r);
    return x != NULL ? x->rate : NAN;
}

const char *pl_network_reaction_source(const pl_network *network, int r) {
    const reaction *x = reaction_at(network, r);
    return x != NULL ? x->source : NULL;
}

int pl_network_reaction_lumped(const pl_network *network, int r) {
    const reaction *x = reaction_at(network, r);
    return x != NULL ? x->lumped : 0;
}

/* ---- the rate law ---- */

/* c^order for a nonnegative integer order, without pow for the common orders. */
static double power(double c, int order) {
    if (order == 0) return 1.0;
    if (order == 1) return c;
    if (order == 2) return c * c;
    return pow(c, order);
}

void pl_network_rates(const pl_network *network, const double *molar, double *rate) {
    for (int r = 0; r < network->reaction_count; r++) {
        const reaction *x = &network->reactions[r];
        double v = x->rate;
        for (int t = 0; t < x->count[LEFT]; t++) v *= power(molar[x->terms[t].species], x->terms[t].coefficient);
        rate[r] = v;
    }
}

/* d[S] += times * (right - left coefficient of S), dynamic species only. */
static void apply(const pl_network *network, const reaction *x, double *d, double times) {
    for (int t = 0; t < x->count[LEFT] + x->count[RIGHT]; t++) {
        int s = x->terms[t].species;
        if (network->species[s].role != PL_DYNAMIC) continue;
        d[s] += (t < x->count[LEFT] ? -times : times) * x->terms[t].coefficient;
    }
}

void pl_network_derivative(const pl_network *network, const double *molar, double *dmolar) {
    for (int s = 0; s < network->species_count; s++) dmolar[s] = 0.0;
    for (int r = 0; r < network->reaction_count; r++) {
        const reaction *x = &network->reactions[r];
        double v = x->rate;
        for (int t = 0; t < x->count[LEFT]; t++) v *= power(molar[x->terms[t].species], x->terms[t].coefficient);
        apply(network, x, dmolar, v);
    }
}

void pl_network_jacobian(const pl_network *network, const double *molar, double *jacobian) {
    int n = network->species_count;
    for (int i = 0; i < n * n; i++) jacobian[i] = 0.0;
    for (int r = 0; r < network->reaction_count; r++) {
        const reaction *x = &network->reactions[r];
        for (int t = 0; t < x->count[LEFT]; t++) {
            int q = x->terms[t].species;
            if (network->species[q].role != PL_DYNAMIC) continue;
            double d = x->rate * x->terms[t].coefficient;  /* d rate / d [q], for this term */
            for (int u = 0; u < x->count[LEFT]; u++)
                d *= power(molar[x->terms[u].species], x->terms[u].coefficient - (u == t));
            for (int v = 0; v < x->count[LEFT] + x->count[RIGHT]; v++) {
                int s = x->terms[v].species;
                if (network->species[s].role != PL_DYNAMIC) continue;
                jacobian[s * n + q] += (v < x->count[LEFT] ? -d : d) * x->terms[v].coefficient;
            }
        }
    }
}
