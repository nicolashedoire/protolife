#ifndef PROTOLIFE_DIAGNOSTIC_H
#define PROTOLIFE_DIAGNOSTIC_H

/*
 * Every protolife call that can fail returns a pl_status. When the caller passes a pl_diagnostic, a failure also
 * says where and why: the code of the library underneath, the model time it happened at, and a message. Nothing is
 * corrected in silence: a call that fails leaves the objects it was given as they were before the call.
 */

typedef enum pl_status {
    PL_OK = 0,
    PL_ERROR_ARGUMENT,     /* an argument is invalid: a null pointer, a size or a tolerance out of range */
    PL_ERROR_MEMORY,       /* an allocation failed */
    PL_ERROR_DOMAIN,       /* a parameter lies outside the range a mechanism was calibrated on */
    PL_ERROR_MODEL,        /* a model callback reported a failure or produced a value that is not finite */
    PL_ERROR_SOLVER,       /* a solver could not reach the requested state */
    PL_ERROR_CONSERVATION, /* a balance drifted beyond its tolerance */
    PL_ERROR_INPUT         /* a file or a text could not be read */
} pl_status;

#define PL_MESSAGE_SIZE 512

typedef struct pl_diagnostic {
    pl_status status;
    int code;                       /* the code of the library underneath (a CVODE flag, ...), 0 when there is none */
    double time;                    /* model time of the failure (s), NAN when it has none */
    char message[PL_MESSAGE_SIZE];  /* empty on success */
} pl_diagnostic;

/* The name of a status, such as "PL_ERROR_SOLVER". */
const char *pl_status_name(pl_status status);

#endif
