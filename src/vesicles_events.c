/* What happens to the vesicles between two integrated stretches (protolife/vesicles.h): ruptures, feeds of micelles,
 * births, divisions, comings apart, the washout and dilutions. After Genesis world/vesicles.c. */
#include <math.h>
#include <stdlib.h>

#include "vesicles_internal.h"

/* A vesicle past its rupture tension lets a share of everything inside out, until its tension is back at the rupture
 * value (found by bisection on the share kept), and reseals. */
void pl_vesicles_rupture(pl_vesicles *p) {
    if (p == NULL) return;
    for (int r = 0; r < p->records; r++) {
        pl_vesicle_shape shape = pl_vesicles_shape(p, r);
        if (!(shape.tension > p->k.rupture)) continue;
        record_t *v = &p->record[r];
        double *own = pl_inside(p, r), osmolarity, osmoles = 0.0;
        pl_vesicles_medium(p, &osmolarity);
        for (int s = 0; s < p->solutes; s++) osmoles += p->solute[s].osmotic * own[s];
        double low = 0.0, high = 1.0;
        for (int i = 0; i < 80; i++) {
            double kept = 0.5 * (low + high);
            if (pl_vesicle_shape_of(&p->k, v->molecules, kept * osmoles, osmolarity).tension > p->k.rupture) high = kept;
            else low = kept;
        }
        for (int s = 0; s < p->solutes; s++) {
            double out = (1.0 - low) * own[s];
            own[s] -= out;
            p->outside[s] += v->count * out;
            p->leaked[s] += v->count * out;
        }
        p->ruptures++;
    }
}

void pl_vesicles_feed_micelles(pl_vesicles *p, double amount) {
    if (p == NULL || !(amount > 0.0)) return;
    const pl_vesicle_constants *k = &p->k;
    double membranes = pl_vesicles_membranes(p), pending = 0.0;
    for (int r = 0; r < p->records; r++) pending += p->record[r].coat;
    double room = fmax(k->coat_share * membranes - pending, 0.0), coated = fmin(amount, room), rest = amount - coated;
    double to_membranes = membranes > 0.0 ? membranes / (membranes + k->nucleation * rest) : 0.0;
    for (int r = 0; r < p->records && membranes > 0.0; r++) {  /* each record its share of the coat and intermediate */
        double share = p->record[r].count * p->record[r].molecules / PL_AVOGADRO / membranes;
        p->record[r].coat += coated * share;
        p->record[r].intermediate += to_membranes * rest * share;
    }
    p->forming += (1.0 - to_membranes) * rest;
    p->nursery = -1;
    p->feeds++;
}

/* What the intermediate made becomes vesicles of the measured radius, closed round the medium (nearly full spheres),
 * joining the record of this feed's newborns; with no record free, it waits. */
void pl_vesicles_births(pl_vesicles *p) {
    const pl_vesicle_constants *k = &p->k;
    if (!(p->newborn > 0.0)) return;
    double molecules = pl_vesicle_molecules_for_radius(k, k->born_radius), count = p->newborn * PL_AVOGADRO / molecules;
    double volume = 0.9 * 4.0 / 3.0 * PL_PI * pow(k->born_radius, 3.0), medium = pl_vesicles_medium(p, NULL);
    double *inside = malloc((size_t)(p->solutes > 0 ? p->solutes : 1) * sizeof *inside);
    if (inside == NULL) return;  /* it waits, as with no record free */
    for (int s = 0; s < p->solutes; s++) inside[s] = 1000.0 * p->outside[s] / medium * volume;
    if (p->nursery < 0 || p->nursery >= p->records) {
        int r = pl_vesicles_append(p, count, molecules, inside, 0);
        if (r < 0) {
            free(inside);
            return;
        }
        p->nursery = r;
        p->record[r].origin = p->feeds;
    } else {  /* more newborns of the same feed: the record keeps their mean, every amount kept */
        record_t *v = &p->record[p->nursery];
        double *own = pl_inside(p, p->nursery), total = v->count + count;
        v->molecules = (v->count * v->molecules + count * molecules) / total;
        v->released = v->count * v->released / total;
        for (int s = 0; s < p->solutes; s++) own[s] = (v->count * own[s] + count * inside[s]) / total;
        v->count = total;
    }
    for (int s = 0; s < p->solutes; s++) p->outside[s] -= count * inside[s];  /* what they enclose leaves the medium */
    free(inside);
    p->born += count;
    p->newborn = 0.0;
}

/* Under enough shear, a vesicle whose reduced volume is below the threshold divides into the fewest equal daughters
 * that reach it, none smaller than the least radius. */
void pl_vesicles_divisions(pl_vesicles *p) {
    const pl_vesicle_constants *k = &p->k;
    if (!(p->shear >= k->critical_shear)) return;
    for (int r = 0; r < p->records; r++) {
        pl_vesicle_shape s = pl_vesicles_shape(p, r);
        double reduced = s.volume / (4.0 / 3.0 * PL_PI * pow(s.radius, 3.0));
        if (!(reduced < k->division_reduced)) continue;
        double n = ceil(pow(k->division_reduced / reduced, 2.0));
        n = fmin(n, floor(pow(s.radius / k->least_radius, 2.0) * (1.0 + 1e-9)));  /* daughters at the least radius */
        if (n < 2.0) continue;
        record_t *v = &p->record[r];
        double *own = pl_inside(p, r);
        p->divided += v->count * (n - 1.0);
        v->count *= n, v->molecules /= n, v->released /= n;
        for (int i = 0; i < p->solutes; i++) own[i] /= n;
        v->generation++;
    }
}

/* Vesicles shrunk below the radius of coming apart: their fatty acid (with what was entering their membranes) joins the
 * monomer, their contents the medium, and their record is taken out. */
void pl_vesicles_come_apart(pl_vesicles *p) {
    for (int r = p->records - 1; r >= 0; r--) {
        record_t *v = &p->record[r];
        double radius = sqrt(v->molecules * p->k.area_per_molecule / 2.0 / (4.0 * PL_PI));
        if (!(radius < p->k.apart_radius)) continue;
        double *own = pl_inside(p, r);
        p->monomer += v->count * v->molecules / PL_AVOGADRO + v->coat + v->intermediate;
        for (int s = 0; s < p->solutes; s++) p->outside[s] += v->count * own[s];
        p->apart += v->count;
        for (int q = r; q + 1 < p->records; q++) {
            p->record[q] = p->record[q + 1];
            double *to = pl_inside(p, q), *from = pl_inside(p, q + 1);
            for (int s = 0; s < p->solutes; s++) to[s] = from[s];
        }
        p->records--;
        if (p->nursery == r) p->nursery = -1;
        else if (p->nursery > r) p->nursery--;
    }
}

/* A share of every record, its micelles and the pools leaves, counted as washed; what left inside vesicles counts as
 * leaving the renewal's ledger. */
static void take_out(pl_vesicles *p, double share) {
    double kept = 1.0 - share;
    for (int r = 0; r < p->records; r++) {
        record_t *v = &p->record[r];
        const double *own = pl_inside(p, r);
        p->washed += share * (v->count * v->molecules / PL_AVOGADRO + v->coat + v->intermediate);
        for (int s = 0; s < p->solutes; s++) p->renewed[s] -= share * v->count * own[s];
        v->count *= kept, v->coat *= kept, v->intermediate *= kept;
    }
    p->washed += share * (p->forming + p->newborn);
    p->forming *= kept, p->newborn *= kept;
}

/* The vesicles that left with the stream over dt, after the step (the medium's renewal is in the step). */
void pl_vesicles_wash(pl_vesicles *p, double dt) {
    if (p->washout > 0.0) take_out(p, 1.0 - exp(-p->washout * dt));
}

void pl_vesicles_dilute(pl_vesicles *p, double share) {
    if (p == NULL || !(share > 0.0) || share > 1.0) return;
    double medium = pl_vesicles_medium(p, NULL);
    take_out(p, share);
    for (int s = 0; s < p->solutes; s++) {
        double fresh = share * p->feed[s] * medium - share * p->outside[s];
        p->outside[s] += fresh, p->renewed[s] += fresh;
    }
    double fresh = share * p->monomer_feed * medium - share * p->monomer;
    p->monomer += fresh, p->renewed_monomer += fresh;
}
