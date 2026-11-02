#include "diagnostic_internal.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>

const char *pl_status_name(pl_status status) {
    switch (status) {
    case PL_OK: return "PL_OK";
    case PL_ERROR_ARGUMENT: return "PL_ERROR_ARGUMENT";
    case PL_ERROR_MEMORY: return "PL_ERROR_MEMORY";
    case PL_ERROR_DOMAIN: return "PL_ERROR_DOMAIN";
    case PL_ERROR_MODEL: return "PL_ERROR_MODEL";
    case PL_ERROR_SOLVER: return "PL_ERROR_SOLVER";
    case PL_ERROR_CONSERVATION: return "PL_ERROR_CONSERVATION";
    case PL_ERROR_INPUT: return "PL_ERROR_INPUT";
    }
    return "PL_UNKNOWN_STATUS";
}

void pl_succeed(pl_diagnostic *diagnostic) {
    if (diagnostic == NULL) return;
    diagnostic->status = PL_OK, diagnostic->code = 0, diagnostic->time = NAN, diagnostic->message[0] = '\0';
}

pl_status pl_fail(pl_diagnostic *diagnostic, pl_status status, int code, double time, const char *format, ...) {
    if (diagnostic == NULL) return status;
    diagnostic->status = status, diagnostic->code = code, diagnostic->time = time;
    va_list arguments;
    va_start(arguments, format);
    vsnprintf(diagnostic->message, sizeof diagnostic->message, format, arguments);
    va_end(arguments);
    return status;
}
