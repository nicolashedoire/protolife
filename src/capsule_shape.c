/* The shape and the elastic energy of a protocell (capsule_internal.h), after Taneja and Higgs's growADE.c. */
#include <math.h>

#include "capsule_internal.h"

#define PI 3.14159265358979323846

capsule_form pl_capsule_form(double volume, double lipid_area) {
    capsule_form f = {0};
    if (volume >= pow(lipid_area, 1.5)) {  /* a swollen sphere: no capsule can hold it */
        f.radius = cbrt(volume), f.length = 0.0, f.actual = pow(volume, 2.0 / 3.0), f.swollen = 1;
        return f;
    }
    /* a capsule of this area and volume: 1.5 r A - 0.5 r^3 = V, by bisection on r as the authors do */
    double low = 0.0, high = sqrt(lipid_area), r;
    do {
        r = 0.5 * (low + high);
        if (1.5 * r * lipid_area - 0.5 * r * r * r >= volume) high = r;
        else low = r;
    } while (high - low > 1e-9);
    f.radius = r, f.length = 2.0 * lipid_area / r - 2.0 * r, f.actual = lipid_area;
    return f;
}

double pl_capsule_energy(const pl_capsule_parameters *k, double DA0, double volume, double outer, double inner) {
    double lipid = 0.5 * (outer + inner);
    capsule_form f = pl_capsule_form(volume, lipid);
    double curvature = k->kappa * PI * (8.0 + f.length / f.radius);
    double shape = f.radius + 0.25 * f.length, held = (outer - inner) / DA0;
    double difference = 16.0 * PI * k->beta * k->kappa / lipid * (shape - held) * (shape - held);
    double actual = fmax(lipid, pow(volume, 2.0 / 3.0)), stretch = k->R0 * (actual - lipid) / k->d;
    double elastic = 4.0 * PI * k->kappa * k->beta * stretch * stretch / lipid;
    return curvature + difference + elastic;
}
