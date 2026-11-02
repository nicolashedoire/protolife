/* Recordings (protolife/scene.h) and the cycles' protocol (protolife/protocol.h). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "protolife/protolife.h"

static char *slurp(const char *path, long *size) {
    FILE *in = fopen(path, "rb");
    if (in == NULL) return NULL;
    fseek(in, 0, SEEK_END);
    *size = ftell(in);
    rewind(in);
    char *text = malloc((size_t)*size + 1);
    if (text != NULL && fread(text, 1, (size_t)*size, in) == (size_t)*size) text[*size] = '\0';
    fclose(in);
    return text;
}

static void test_scene(void) {
    pl_scene *s = NULL;
    pl_diagnostic d;
    CHECK(pl_scene_create("lineage", "A \"quoted\" title", "", "tests", &s, &d) == PL_OK);
    CHECK(pl_scene_add_value(s, "radius", "radius", "nm", &d) == PL_OK);
    CHECK(pl_scene_add_value(s, "count", "count", "", &d) == PL_OK);
    CHECK(pl_scene_add_constant(s, "cycles", "cycles", "", 4.0, &d) == PL_OK);
    CHECK(pl_scene_open(s, PROTOLIFE_TEST_DIR "/scene-test/nested", &d) == PL_OK);
    CHECK(pl_scene_add_value(s, "late", "late", "", &d) == PL_ERROR_ARGUMENT);  /* the manifest is fixed */
    float row[2] = {500.0f, 1.0f};
    for (int k = 0; k < 3; k++) CHECK(pl_scene_frame(s, 10.0 * k, row, NULL, &d) == PL_OK);
    CHECK(pl_scene_frames(s) == 3 && pl_scene_value_count(s) == 2);
    pl_scene_free(s);
    long size = 0;
    char *json = slurp(PROTOLIFE_TEST_DIR "/scene-test/nested/scene.json", &size);
    CHECK(json != NULL && strstr(json, "\"format\": \"protolife-scene-1\"") != NULL);
    CHECK(json != NULL && strstr(json, "A \\\"quoted\\\" title") != NULL && strstr(json, "\"frame_bytes\": 16") != NULL);
    free(json);
    char *frames = slurp(PROTOLIFE_TEST_DIR "/scene-test/nested/frames.bin", &size);
    CHECK(frames != NULL && size == 3 * 16);
    free(frames);
}

static void test_cycles(void) {
    pl_vesicle_constants k = pl_vesicle_oleate();
    pl_vesicles *p = NULL;
    CHECK(pl_vesicles_create(&k, 0, &p, NULL) == PL_OK);
    double molecules = pl_vesicle_molecules_for_radius(&k, 100e-9);
    CHECK(pl_vesicles_add(p, 1e15, molecules, NULL, 1, NULL, NULL) == PL_OK);
    pl_cycles c = {.count = 2, .start = 100.0, .period = 60.0, .dilution = 0.5, .micelles = 1.0, .feeds = 1,
                   .interval = 0.0, .shear = 20.0, .stirred = 0};
    double membranes = pl_vesicles_membranes(p);
    pl_cycles_apply(p, &c, 99.0, 100.0);  /* the first cycle starts: half the vessel replaced */
    CHECK_CLOSE(pl_vesicles_record(p, 0).count, 0.5e15, 1e-12);
    pl_cycles_apply(p, &c, 109.0, 110.0);  /* 10 s in: one equivalent of micelles */
    CHECK_CLOSE(pl_vesicles_entering(p), 0.5 * membranes, 1e-12);
    pl_cycles_apply(p, &c, 139.0, 140.0);  /* two thirds in: shear for 30 s */
    CHECK(pl_vesicles_shear(p) == 20.0);
    pl_cycles_apply(p, &c, 169.0, 170.0);
    CHECK(pl_vesicles_shear(p) == 0.0);
    CHECK_CLOSE(pl_cycles_kept(&c, 165.0), 0.25, 1e-15);  /* two dilutions by then */
    pl_vesicles_free(p);
}

int main(void) {
    test_scene();
    test_cycles();
    return check_report("scene");
}
