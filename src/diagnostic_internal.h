#ifndef PROTOLIFE_DIAGNOSTIC_INTERNAL_H
#define PROTOLIFE_DIAGNOSTIC_INTERNAL_H

#include "protolife/diagnostic.h"

#if defined(__GNUC__) || defined(__clang__)
#define PL_PRINTF(pattern, first) __attribute__((format(printf, pattern, first)))
#else
#define PL_PRINTF(pattern, first)
#endif

/* Marks a success: status PL_OK, no code, no time, no message. A null diagnostic is allowed. */
void pl_succeed(pl_diagnostic *diagnostic);

/* Fills the diagnostic, when there is one, and returns `status`. */
pl_status pl_fail(pl_diagnostic *diagnostic, pl_status status, int code, double time, const char *format, ...)
    PL_PRINTF(5, 6);

#endif
