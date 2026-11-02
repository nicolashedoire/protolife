/* The reversible Oregonator (Belousov-Zhabotinsky) in a closed vessel, printed as CSV: time (s), then the
 * concentration (mol/L) of every dynamic species, and the integrator's cumulative step count.
 *
 *   oregonator [network.net] [duration_s] [interval_s]
 */
#include <stdio.h>
#include <stdlib.h>

#include "protolife/protolife.h"

int main(int argc, char **argv) {
    const char *path = argc > 1 ? argv[1] : PROTOLIFE_DATA "/networks/oregonator-reversible.net";
    double duration = argc > 2 ? atof(argv[2]) : 60.0, interval = argc > 3 ? atof(argv[3]) : 0.05;
    pl_network *network = NULL;
    pl_reactor *reactor = NULL;
    pl_diagnostic d;
    if (pl_network_load(path, &network, &d) != PL_OK || pl_network_check_balance(network, &d) != PL_OK) {
        fprintf(stderr, "%s: %s\n", pl_status_name(d.status), d.message);
        pl_network_free(network);
        return 1;
    }
    pl_reactor_options options = {.relative_tolerance = 1e-11, .absolute_tolerance = 1e-18, .max_steps = 1000000};
    if (pl_reactor_create(network, &options, &reactor, &d) != PL_OK) {
        fprintf(stderr, "%s: %s\n", pl_status_name(d.status), d.message);
        pl_network_free(network);
        return 1;
    }
    int n = pl_network_species_count(network);
    printf("time");
    for (int s = 0; s < n; s++)
        if (pl_network_species_role(network, s) == PL_DYNAMIC) printf(",%s", pl_network_species_name(network, s));
    printf(",steps\n");
    int status = 0;
    for (long k = 0; status == 0 && k * interval <= duration + 1e-9; k++) {
        if (k > 0 && pl_reactor_advance(reactor, k * interval, &d) != PL_OK) {
            fprintf(stderr, "%s: %s\n", pl_status_name(d.status), d.message);
            status = 1;
            break;
        }
        printf("%.6g", k * interval);
        for (int s = 0; s < n; s++)
            if (pl_network_species_role(network, s) == PL_DYNAMIC) printf(",%.9g", pl_reactor_molar(reactor)[s]);
        pl_ode_statistics statistics;
        pl_reactor_get_statistics(reactor, &statistics);
        printf(",%ld\n", statistics.steps);
    }
    pl_reactor_free(reactor);
    pl_network_free(network);
    return status;
}
