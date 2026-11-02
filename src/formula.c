/* Chemical formulas (formula_internal.h), after Genesis's chemistry/composition.c, extended to every element. */
#include "formula_internal.h"

#include <ctype.h>
#include <string.h>

static const char *const SYMBOLS[PL_ELEMENTS] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
    "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

const char *pl_element_symbol(int element) {
    return element >= 0 && element < PL_ELEMENTS ? SYMBOLS[element] : "?";
}

static int element_find(const char *symbol, size_t length) {
    for (int e = 0; e < PL_ELEMENTS; e++)
        if (strlen(SYMBOLS[e]) == length && strncmp(SYMBOLS[e], symbol, length) == 0) return e;
    return -1;
}

/* A count after an element or a group: 1 when absent; -1 when it overflows. */
static int read_count(const char **p) {
    if (!isdigit((unsigned char)**p)) return 1;
    long n = 0;
    while (isdigit((unsigned char)**p)) {
        n = 10 * n + (*(*p)++ - '0');
        if (n > 1000000) return -1;
    }
    return n > 0 ? (int)n : -1;
}

/* Elements and parenthesised groups until ')' or the end, added `times` times; -1 on error. */
static int parse_group(const char **p, int atoms[PL_ELEMENTS], int depth) {
    while (**p != '\0' && **p != ')') {
        int part[PL_ELEMENTS] = {0};
        if (**p == '(') {
            (*p)++;
            if (depth > 8 || parse_group(p, part, depth + 1) != 0 || **p != ')') return -1;
            (*p)++;
        } else if (isupper((unsigned char)**p)) {
            size_t length = islower((unsigned char)(*p)[1]) ? 2 : 1;
            int element = element_find(*p, length);
            if (element < 0) return -1;
            part[element] = 1;
            *p += length;
        } else {
            return -1;
        }
        int times = read_count(p);
        if (times < 0) return -1;
        for (int e = 0; e < PL_ELEMENTS; e++) atoms[e] += times * part[e];
    }
    return 0;
}

int pl_formula_parse(const char *formula, int atoms[PL_ELEMENTS]) {
    memset(atoms, 0, PL_ELEMENTS * sizeof *atoms);
    if (formula == NULL || *formula == '\0') return -1;
    const char *p = formula;
    return parse_group(&p, atoms, 0) == 0 && *p == '\0' ? 0 : -1;
}
