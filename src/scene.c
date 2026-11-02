/* Recordings of simulations for a viewer (protolife/scene.h), after Genesis world/scene_writer.c. */
#include "protolife/scene.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#define make_directory(path) _mkdir(path)
#else
#define make_directory(path) mkdir(path, 0755)
#endif

#include "diagnostic_internal.h"

#define SCENE_FORMAT "protolife-scene-1"

typedef struct {
    char *key, *label, *unit;
    double value;
} quantity;

typedef struct {
    quantity *item;
    int count, capacity;
} quantities;

struct pl_scene {
    char *kind, *title, *subtitle, *source, *time_unit;
    int n;
    double box;
    quantities values, fields, constants;
    FILE *frames;
    long cells, written;
};

static char *copy_text(const char *text) {
    if (text == NULL) text = "";
    size_t n = strlen(text) + 1;
    char *copy = malloc(n);
    if (copy != NULL) memcpy(copy, text, n);
    return copy;
}

static void free_quantities(quantities *q) {
    for (int i = 0; i < q->count; i++) free(q->item[i].key), free(q->item[i].label), free(q->item[i].unit);
    free(q->item);
}

void pl_scene_free(pl_scene *s) {
    if (s == NULL) return;
    if (s->frames != NULL) fclose(s->frames);
    free(s->kind), free(s->title), free(s->subtitle), free(s->source), free(s->time_unit);
    free_quantities(&s->values), free_quantities(&s->fields), free_quantities(&s->constants);
    free(s);
}

pl_status pl_scene_create(const char *kind, const char *title, const char *subtitle, const char *source,
                          pl_scene **out, pl_diagnostic *diagnostic) {
    if (out == NULL || kind == NULL || *kind == '\0')
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "a scene needs a kind and an address");
    *out = NULL;
    pl_scene *s = calloc(1, sizeof *s);
    if (s == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "the scene could not be allocated");
    s->n = 1;
    s->kind = copy_text(kind), s->title = copy_text(title), s->subtitle = copy_text(subtitle), s->source = copy_text(source);
    s->time_unit = copy_text("s");
    if (s->kind == NULL || s->title == NULL || s->subtitle == NULL || s->source == NULL || s->time_unit == NULL) {
        pl_scene_free(s);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "the scene could not be allocated");
    }
    *out = s;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_scene_set_grid(pl_scene *s, int n, double box, pl_diagnostic *diagnostic) {
    if (s == NULL || s->frames != NULL || n < 1 || !(box >= 0.0) || !isfinite(box))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "a grid of at least one cell, before the scene is opened");
    s->n = n, s->box = box;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_scene_set_time_unit(pl_scene *s, const char *unit, pl_diagnostic *diagnostic) {
    if (s == NULL || unit == NULL || *unit == '\0' || s->frames != NULL)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "a time unit, before the scene is opened");
    char *copy = copy_text(unit);
    if (copy == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory for the time unit");
    free(s->time_unit);
    s->time_unit = copy;
    pl_succeed(diagnostic);
    return PL_OK;
}

static pl_status add(pl_scene *s, quantities *q, const char *key, const char *label, const char *unit, double value,
                     pl_diagnostic *diagnostic) {
    if (s == NULL || key == NULL || *key == '\0' || s->frames != NULL)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "a quantity needs a key, before the scene is opened");
    if (q->count == q->capacity) {
        int capacity = q->capacity > 0 ? 2 * q->capacity : 16;
        quantity *bigger = realloc(q->item, (size_t)capacity * sizeof *bigger);
        if (bigger == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory for the quantity %s", key);
        q->item = bigger, q->capacity = capacity;
    }
    quantity x = {copy_text(key), copy_text(label), copy_text(unit), value};
    if (x.key == NULL || x.label == NULL || x.unit == NULL) {
        free(x.key), free(x.label), free(x.unit);
        return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory for the quantity %s", key);
    }
    q->item[q->count++] = x;
    pl_succeed(diagnostic);
    return PL_OK;
}

pl_status pl_scene_add_value(pl_scene *s, const char *key, const char *label, const char *unit, pl_diagnostic *d) {
    return add(s, s != NULL ? &s->values : NULL, key, label, unit, 0.0, d);
}

pl_status pl_scene_add_field(pl_scene *s, const char *key, const char *label, const char *unit, pl_diagnostic *d) {
    return add(s, s != NULL ? &s->fields : NULL, key, label, unit, 0.0, d);
}

pl_status pl_scene_add_constant(pl_scene *s, const char *key, const char *label, const char *unit, double value,
                                pl_diagnostic *d) {
    if (!isfinite(value)) return pl_fail(d, PL_ERROR_ARGUMENT, 0, NAN, "the constant %s is not finite", key != NULL ? key : "");
    return add(s, s != NULL ? &s->constants : NULL, key, label, unit, value, d);
}

static int make_directories(const char *path) {
    size_t n = strlen(path);
    char *partial = malloc(n + 1);
    if (partial == NULL) return -1;
    memcpy(partial, path, n + 1);
    int ok = 1;
    for (char *p = partial + 1; ok && *p != '\0'; p++) {
        if (*p != '/' && *p != '\\') continue;
        char separator = *p;
        *p = '\0';
        ok = make_directory(partial) == 0 || errno == EEXIST;
        *p = separator;
    }
    ok = ok && (make_directory(partial) == 0 || errno == EEXIST);
    free(partial);
    return ok ? 0 : -1;
}

static void json_string(FILE *out, const char *s) {
    fputc('"', out);
    for (const unsigned char *c = (const unsigned char *)s; *c != '\0'; c++) {
        if (*c == '"' || *c == '\\') fprintf(out, "\\%c", *c);
        else if (*c == '\n') fputs("\\n", out);
        else if (*c < 0x20) fprintf(out, "\\u%04x", *c);
        else fputc(*c, out);
    }
    fputc('"', out);
}

static void write_quantities(FILE *out, const char *name, const quantities *q, int with_values) {
    fprintf(out, "  \"%s\": [", name);
    for (int i = 0; i < q->count; i++) {
        fputs(i ? ",\n    {\"key\": " : "\n    {\"key\": ", out);
        json_string(out, q->item[i].key);
        fputs(", \"label\": ", out);
        json_string(out, q->item[i].label);
        fputs(", \"unit\": ", out);
        json_string(out, q->item[i].unit);
        if (with_values) fprintf(out, ", \"value\": %.17g", q->item[i].value);
        fputc('}', out);
    }
    fputs(q->count ? "\n  ]" : "]", out);
}

pl_status pl_scene_open(pl_scene *s, const char *directory, pl_diagnostic *diagnostic) {
    if (s == NULL || directory == NULL || s->frames != NULL)
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, NAN, "a scene not yet opened and a directory are required");
    if (make_directories(directory) != 0) return pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "cannot make %s", directory);
    size_t length = strlen(directory) + 16;
    char *path = malloc(length);
    if (path == NULL) return pl_fail(diagnostic, PL_ERROR_MEMORY, 0, NAN, "no memory to open the scene");
    s->cells = (long)s->n * s->n;
    long bytes = (long)sizeof(double) + (long)sizeof(float) * (s->values.count + s->fields.count * s->cells);
    snprintf(path, length, "%s/scene.json", directory);
    FILE *out = fopen(path, "w");
    if (out == NULL) {
        pl_status status = pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "cannot write %s", path);
        free(path);
        return status;
    }
    fprintf(out, "{\n  \"format\": \"%s\",\n  \"kind\": ", SCENE_FORMAT);
    json_string(out, s->kind);
    fputs(",\n  \"title\": ", out);
    json_string(out, s->title);
    fputs(",\n  \"subtitle\": ", out);
    json_string(out, s->subtitle);
    fputs(",\n  \"source\": ", out);
    json_string(out, s->source);
    fprintf(out, ",\n  \"grid\": {\"n\": %d, \"box_m\": %.17g},\n  \"time_unit\": ", s->n, s->box);
    json_string(out, s->time_unit);
    fputs(",\n", out);
    fprintf(out, "  \"frame_bytes\": %ld,\n  \"frame_layout\": \"float64 time, float32 values, float32 fields x cells\",\n",
            bytes);
    write_quantities(out, "values", &s->values, 0);
    fputs(",\n", out);
    write_quantities(out, "fields", &s->fields, 0);
    fputs(",\n", out);
    write_quantities(out, "constants", &s->constants, 1);
    fputs("\n}\n", out);
    int written = fclose(out) == 0;
    snprintf(path, length, "%s/frames.bin", directory);
    s->frames = written ? fopen(path, "wb") : NULL;
    pl_status status = s->frames != NULL ? PL_OK : pl_fail(diagnostic, PL_ERROR_INPUT, 0, NAN, "cannot write %s", path);
    free(path);
    if (status == PL_OK) pl_succeed(diagnostic);
    return status;
}

pl_status pl_scene_frame(pl_scene *s, double time, const float *values, const float *fields, pl_diagnostic *diagnostic) {
    if (s == NULL || s->frames == NULL) return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, time, "the scene is not open");
    if ((s->values.count > 0 && values == NULL) || (s->fields.count > 0 && fields == NULL))
        return pl_fail(diagnostic, PL_ERROR_ARGUMENT, 0, time, "a frame needs its values and fields");
    size_t field_floats = (size_t)s->fields.count * (size_t)s->cells, value_floats = (size_t)s->values.count;
    int ok = fwrite(&time, sizeof time, 1, s->frames) == 1 &&
             (value_floats == 0 || fwrite(values, sizeof *values, value_floats, s->frames) == value_floats) &&
             (field_floats == 0 || fwrite(fields, sizeof *fields, field_floats, s->frames) == field_floats) &&
             fflush(s->frames) == 0;
    if (!ok) return pl_fail(diagnostic, PL_ERROR_INPUT, 0, time, "could not append a frame");
    s->written++;
    pl_succeed(diagnostic);
    return PL_OK;
}

int pl_scene_value_count(const pl_scene *s) { return s != NULL ? s->values.count : 0; }
long pl_scene_frames(const pl_scene *s) { return s != NULL ? s->written : 0; }
