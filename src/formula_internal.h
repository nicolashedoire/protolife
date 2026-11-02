#ifndef PROTOLIFE_FORMULA_INTERNAL_H
#define PROTOLIFE_FORMULA_INTERNAL_H

/* Atoms of each element in a chemical formula, over the whole periodic table. */
#define PL_ELEMENTS 118

/* The symbol of element e (0 is hydrogen). */
const char *pl_element_symbol(int element);

/* The index of an element from its symbol, or -1. */
int pl_element_find(const char *symbol);

/* Reads a formula such as "HBrO2" or "CH2(COOH)2" into atoms[PL_ELEMENTS]: 0, or -1 with `atoms` unspecified. */
int pl_formula_parse(const char *formula, int atoms[PL_ELEMENTS]);

#endif
