/* Uses an installed protolife: loads a network it ships and runs Taneja and Higgs's cell for a while. */
#include <stdio.h>

#include "protolife/protolife.h"

int main(void) {
    pl_network *network = NULL;
    pl_capsule *cell = NULL;
    pl_diagnostic d;
    double mu[3] = {0.01, 0.0, 0.015}, outside[3] = {1.0, 0.0, 0.0};
    pl_capsule_parameters k = pl_capsule_taneja_higgs();
    if (pl_network_load(NETWORKS "/taneja-higgs.net", &network, &d) != PL_OK ||
        pl_capsule_create(&k, network, mu, outside, 1, &cell, &d) != PL_OK ||
        pl_capsule_advance(cell, 100.0, 0.01, &d) != PL_OK) {
        fprintf(stderr, "%s: %s\n", pl_status_name(d.status), d.message);
        return 1;
    }
    printf("protolife %s works: after 100 time units the cell holds %.4g V0\n", PL_VERSION, pl_capsule_get_shape(cell).volume);
    pl_capsule_free(cell);
    pl_network_free(network);
    return 0;
}
