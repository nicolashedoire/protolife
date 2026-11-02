# Architecture

This is the design the code will follow. It is a proposal until the C core exists.

## Layers

| Layer | Contents | Origin |
|---|---|---|
| **numerics** | stiff integration with events, linear solvers, FFT, random numbers, stochastic simulation | SUNDIALS and pocketfft as dependencies; Gillespie and hybrid simulation from Genesis |
| **chemistry** | molecules as graphs (SMILES in and out), stereochemistry, reaction rules, network growth selected by flux, thermochemistry by group additivity, kinetics from barriers, autocatalytic cores and chemical organisations | Genesis, with RMG's core/edge criterion and RDKit's canonical ranking borrowed |
| **compartments** | vesicles (osmosis, tension, rupture, micelle feeding, shear division), droplets and condensates | Genesis, checked against Taneja & Higgs 2025 and Mavelli & Ruiz-Mirazo |
| **media and protocols** | a well-stirred vessel, a flow reactor, a dialysis bag, grids for reaction–diffusion; feeding, dilution and shear cycles | Genesis |
| **observation** | matter balances, lineages, signatures of a chemistry, controls (lifeless twin, random parents) | Genesis, with StochPy's lineage and population views borrowed |
| **exchange** | SBML Level 3 Version 2 export (mass-action subset, SMILES and InChI as annotations), Cantera YAML, the CatReNet text format, a scene format for viewers | new |
| **python** | `import protolife`: the objects above, NumPy arrays, typed exceptions, plotting helpers | new |

## C core

- **C11**, built with CMake as static and shared libraries.
- **Opaque handles** and a stable C API under one prefix (`pl_`). The Python layer, and
  any other language, binds to that API only.
- **SI units** everywhere in the API; each quantity documents its unit.
- **Struct of arrays**: each per-vesicle or per-species field is one contiguous
  `double` array, so Python sees it without a copy.
- **Status codes**: every call that can fail returns a code and fills a diagnostic
  (time, variable, reason). Nothing is clamped or corrected in silence.
- **Validity domains**: each mechanism carries its source and the parameter ranges
  where it was calibrated; a value outside them is an error unless the caller
  explicitly declares an extrapolation.

## Python layer

- **nanobind** bindings over the C API, built by **scikit-build-core**, wheels for
  Linux, macOS and Windows by **cibuildwheel**, published on PyPI then conda-forge.
- Arrays are exchanged as NumPy arrays without copies. A wrong dtype or layout is
  rejected, never silently converted.
- Each C status maps to a Python exception, for example
  `SolverFailure(code, time, diagnostics)`.
- Optional extras: `protolife[rdkit]` for drawings, reaction SMARTS and RDKit
  conversion; `protolife[sbml]` for libSBML validation.

## Quality

- C unit tests with CTest, under AddressSanitizer and UndefinedBehaviorSanitizer.
- Python tests with pytest, run on every wheel.
- **Regression**: each brick taken from Genesis reproduces the Genesis result it comes
  from, within a stated tolerance (exactly where the numerics are unchanged).
- **Reference checks** in separate CI jobs: libRoadRunner for exported networks,
  GillesPy2 for stochastic runs, LAMMPS for coarse-grained lipids, MØD for formose
  networks, CatReNet for autocatalytic sets. These tools are never shipped.
- **Published models as tests**: Taneja & Higgs 2025 (growth–division cycle only when
  lipid supply exceeds equilibrium; division lost above a flip-flop rate), Mavelli &
  Ruiz-Mirazo's stationary size, Villani & Serra's synchronisation.

## Documentation

Sphinx with the PyData theme, numpydoc docstrings, and a sphinx-gallery of runnable
examples that double as notebooks. Each mechanism page states its equations, its
source, its validity domain and the test that checks it.

## Integration steps

1. Landscape of existing tools (done).
2. C core: numerics on SUNDIALS, then chemistry, compartments, media, observation,
   one brick at a time, each with its regression test.
3. Python layer and the first example (two lineages, one with chemistry, one without).
4. Documentation, catalogue of calibrated worlds, reference tests against published
   models, first release on PyPI.
