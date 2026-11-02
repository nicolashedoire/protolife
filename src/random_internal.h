#ifndef PROTOLIFE_RANDOM_INTERNAL_H
#define PROTOLIFE_RANDOM_INTERNAL_H

#include <math.h>
#include <stdint.h>

/* splitmix64, then 53 bits: a uniform number in [0, 1), the same on every platform. */
static inline double pl_uniform(uint64_t *state) {
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return (double)((z ^ (z >> 31)) >> 11) * 0x1.0p-53;
}

/* A draw of Binomial(n, p): n Bernoulli trials up to 1000, beyond that the normal approximation, rounded and kept
 * within [0, n] (its error is below the noise of the draw itself at such numbers). */
static inline long long pl_binomial(uint64_t *state, long long n, double p) {
    if (n <= 0 || !(p > 0.0)) return 0;
    if (p >= 1.0) return n;
    if (n <= 1000) {
        long long k = 0;
        for (long long i = 0; i < n; i++) k += pl_uniform(state) < p;
        return k;
    }
    double u1 = pl_uniform(state), u2 = pl_uniform(state);
    double normal = sqrt(-2.0 * log(1.0 - u1)) * cos(2.0 * 3.14159265358979323846 * u2);
    double k = (double)llround((double)n * p + normal * sqrt((double)n * p * (1.0 - p)));
    return k < 0 ? 0 : k > (double)n ? n : (long long)k;
}

#endif
