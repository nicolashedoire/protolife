#ifndef PROTOLIFE_CHECK_H
#define PROTOLIFE_CHECK_H

/* A minimal test harness: each test program counts its failed checks and exits nonzero when there is one. */
#include <math.h>
#include <stdio.h>

static int check_failures;

static inline void check_true(int ok, const char *what, const char *file, int line) {
    if (ok) return;
    fprintf(stderr, "%s:%d: failed: %s\n", file, line, what);
    check_failures++;
}

/* |actual - expected| <= relative * |expected| */
static inline void check_close(double actual, double expected, double relative, const char *what, const char *file,
                        int line) {
    if (fabs(actual - expected) <= relative * fabs(expected)) return;
    fprintf(stderr, "%s:%d: failed: %s: %.17g against %.17g (relative %.3g, allowed %.3g)\n", file, line, what, actual,
            expected, fabs(actual - expected) / fabs(expected), relative);
    check_failures++;
}

#define CHECK(condition) check_true((condition), #condition, __FILE__, __LINE__)
#define CHECK_CLOSE(actual, expected, relative) \
    check_close((actual), (expected), (relative), #actual, __FILE__, __LINE__)

static inline int check_report(const char *name) {
    if (check_failures == 0) printf("%s: all checks passed\n", name);
    else fprintf(stderr, "%s: %d failed checks\n", name, check_failures);
    return check_failures == 0 ? 0 : 1;
}

#endif
