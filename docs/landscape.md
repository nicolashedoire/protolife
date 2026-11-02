# Landscape of existing tools

Checked on 2 November 2026 from the projects' own pages (repositories, PyPI, papers).
Facts we could not confirm are marked *unverified*. Each tool gets one verdict:

- **depend**: protolife uses it, as a required or optional dependency;
- **interoperate**: protolife exchanges files with it, or checks its results against it;
- **borrow**: protolife re-implements an equation, algorithm, data set or API idea from it,
  with credit;
- **ignore**: not useful to protolife, for the reason given.

protolife is BSD-3-Clause, so GPL code can only sit *next to* it (separate processes,
file exchange, reference checks in CI), never under it.

## Conclusion

No maintained, permissively licensed library simulates a population of vesicles with
physical membrane mechanics and internal chemistry. Existing protocell code is either a
teaching toy or tied to a single paper. The closest published models are worth
re-implementing as reference tests. For numerics, chemistry identifiers and packaging,
mature permissive tools exist and protolife should use them instead of hand-written code.

## Protocell and origin-of-life simulators

| Tool | What it simulates | Language | License | Last activity | Verdict |
|---|---|---|---|---|---|
| [Taneja & Higgs 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12113545/), [`growADE.c`](https://zenodo.org/records/15157714) | Osmotic growth dV/dt = λ·A·ΔC, two leaflets with lipid exchange and flip-flop, division when the area-difference elastic energy drops, bursting, an autocatalytic reaction inside | C | CC-BY-4.0 | 2025-04 | **borrow**: reference case and regression targets; re-implement rather than copy (CC-BY is not a code license) |
| ENVIRONMENT (Mavelli & Ruiz-Mirazo, [2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC3303737/)–[2013](https://pubs.rsc.org/en/Content/ArticleLanding/2013/IB/C2IB20222K)) | Gillespie with a changing volume; membrane area from lipid count; burst and division bounds on the shape factor; lipid uptake and release; the osmotic synchronisation condition | — | no public code found | papers | **borrow**: stability window, lipid kinetics, analytic stationary size as a test |
| Villani, Serra et al., [Life 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4284470/) | Container growth coupled to a catalytic network of binary strings; division into halves at a threshold | — | no code released (*unverified*) | paper | **borrow**: "generation times synchronise only if an autocatalytic set is coupled to growth" as a property test |
| [StochPy](https://github.com/SystemsBioinformatics/stochpy) | Gillespie and tau-leap, SBML; a cell-division module with binomial partitioning and lineage vs extant-population views | Python | BSD-3-Clause | 2025-06 | **borrow**: the division and lineage API |
| [ProtoSim](https://github.com/Piltxi/ProtoSim) | ODE protocell, membrane amount plus internal species, division at a threshold; no osmosis, no units | Python | none stated | 2024-05 | **ignore**: no license; the model comes from Villani & Serra |
| [lifesimulatoR](https://noushinn.github.io/lifesimulatoR/) | Abundance grows, leaks, divides above a threshold; no membrane | R | MIT | 2026-06 | **ignore**: teaching toy |
| [ARCHNET](https://github.com/segrelab/string-chemistry) ([paper](https://link.springer.com/article/10.1007/s00239-021-10018-0)) | String chemistries turned into flux-balance models | Python | none stated | 2026-01 | **interoperate**: compare flux-based selection via SBML |
| [GARD](https://github.com/ModelingOriginsofLife/GARD) | Lipid-world composomes, compositional heredity | MATLAB | none stated | 2014 | **borrow**: the composome similarity idea (exact formula *unverified*) |
| [PyCellChemistry](https://www.cs.mun.ca/~banzhaf/AC-home.html) | Nested compartments under a hierarchical Gillespie | Python | GPL-3.0 | *unverified* | **borrow**: the nested-compartment API idea only |
| [PyAChemKit](https://github.com/afaulconbridge/PyAChemKit), [Squirm3](https://github.com/timhutton/squirm3), [Organic Builder](https://sourceforge.net/projects/organicbuilder/) | Abstract artificial chemistries | Python 2 / C++ / Java | BSD / GPL-3.0 / GPL-2.0 | 2011–2018 | **ignore**: dormant, abstract |

## Chemistry, reaction networks and thermodynamics

| Tool | What matters to protolife | Language | License | Last release | Verdict |
|---|---|---|---|---|---|
| [RDKit](https://github.com/rdkit/rdkit) | SMILES/SMARTS, full stereochemistry, reaction SMARTS with chirality handling | C++ / Python | BSD-3-Clause | 2026.09.1 | **depend (optional Python extra)** and **borrow** its stereo-aware canonical ranking ([Schneider, Sayle & Landrum 2015](https://doi.org/10.1021/acs.jcim.5b00543)); test oracle for our canonical SMILES |
| [InChI](https://github.com/IUPAC-InChI/InChI) | Canonical identifier with stereo layers, plain C | C | MIT | 1.07.5 | **depend (optional, C core)** for export as an exchange key |
| [Indigo](https://github.com/epam/Indigo) | Plain C API, canonical SMILES, stereo, reaction enumeration | C/C++ | Apache-2.0 | 1.46.0 | fallback if a C-callable stereo backend is ever needed |
| [MØD](https://jakobandersen.github.io/mod/) | Chemistry as graph rewriting, derivation graphs, hyperflows for autocatalysis | C++ / Python | GPL-3.0 | 1.0.0 (2025-05) | **interoperate** (benchmark: build the same formose network in both) and **borrow** rule semantics; never link |
| [RMG-Py](https://github.com/ReactionMechanismGenerator/RMG-Py) | Rule-based mechanism growth selected by flux, group-additivity thermo, rate rules | Python | MIT | 4.0.0 | **borrow**: the core/edge flux criterion and hierarchical group trees |
| [Cantera](https://cantera.org/) | Thermo and kinetics, YAML mechanisms with inline units, generated C API | C++ / Python | BSD-3-Clause | 3.2.0 | **borrow** the API design; **interoperate** through YAML export |
| [libSBML](https://github.com/sbmlteam/libsbml) | SBML read/write | C++ | LGPL-2.1-or-later | 5.21.2 | optional Python dependency only, or write the XML ourselves |
| [libRoadRunner](https://pypi.org/project/libroadrunner/), [Tellurium](https://pypi.org/project/tellurium/), [COPASI](https://github.com/copasi/COPASI) | SBML simulation and analysis | C++ / Python | Apache-2.0 / Apache-2.0 / Artistic-2.0 | 2026 | **interoperate**: libRoadRunner as a reference simulator in CI |
| [CatReNet](https://github.com/husonlab/catrenet) | Maximal RAF, CAF and pseudo-RAF sets | Java | GPL-3.0 (one source reads "other", *unverified*) | 1.1.4 | **interoperate**: export its text format to cross-check autocatalytic cores |
| [PySB](https://github.com/pysb/pysb) / [BioNetGen](https://github.com/RuleWorld/bionetgen) | Rule-based models on site-graphs | Python / C++ | BSD-2-Clause / MIT | 2026 / 2025 | **borrow** PySB's model-building style |
| [eQuilibrator](https://equilibrator.readthedocs.io/) | Aqueous ΔG'° with uncertainty, pH and ionic-strength corrections | Python | MIT (data license *unverified*) | 2026 | **borrow**: check our group additivity on sugars it knows |
| [Open Babel](https://github.com/openbabel/openbabel) | Format conversion | C++ | GPL-2.0 | 3.2.1 | **ignore**: GPL, RDKit covers it |

## Numerics and other scales

| Tool | What it gives | Language | License | Last release | Verdict |
|---|---|---|---|---|---|
| [SUNDIALS](https://github.com/llnl/sundials) | Stiff integrators (CVODE, ARKODE), GMRES, root-finding for events, error codes returned to the caller | C | BSD-3-Clause | 7.9.0 | **depend**: replaces hand-written integrators; build without KLU (LGPL) |
| [pocketfft](https://github.com/mreineck/pocketfft) / [KissFFT](https://github.com/mborgerding/kissfft) | FFT | C | BSD-3-Clause | — / 131.2.0 | **depend** (vendored); never FFTW, which is GPL |
| [GillesPy2](https://pypi.org/project/gillespy2/) | Stochastic simulation | Python / C++ | GPL-3.0 | 1.8.3 | **interoperate**: reference for our Gillespie |
| [LAMMPS](https://docs.lammps.org/pair_cosine_squared.html) | Ships the Cooke lipid potential | C++ | GPL-2.0 | 2025 | **interoperate**: reference for coarse-grained lipids |
| [Smoldyn](https://github.com/ssandrews/Smoldyn), [ReaDDy](https://github.com/readdy/readdy), [MCell](https://mcell.org/) | Particle reaction–diffusion near surfaces | C/C++ | LGPL or GPL (sources disagree) / BSD-3-Clause / MIT | 2026 | **interoperate** (references); **borrow** ReaDDy's topology design |
| [HOOMD-blue](https://github.com/glotzerlab/hoomd-blue), [OpenMM](https://github.com/openmm/openmm) | GPU molecular dynamics | C++ | BSD-3-Clause / MIT + LGPL | 2026 | **interoperate**: possible GPU route for large membrane runs |
| [Morpheus](https://morpheus.gitlab.io/), [Chaste](https://github.com/Chaste/Chaste) | Cells with internal ODEs and division | C++ | BSD-3-Clause | 2025 / 2026 | **borrow**: declarative models (Morpheus), per-cell model plug-ins (Chaste) |
| [PhysiCell](https://github.com/MathCancer/PhysiCell), [CompuCell3D](https://github.com/CompuCell3D/CompuCell3D) | Tissue-scale agents | C++ | BSD-3-Clause / MIT | 2025 | **ignore**: tissue scale, not protocell chemistry |

## Packaging: how C cores reach Python users

- [SUNDIALS' own Python bindings](https://sundials.readthedocs.io/en/latest/developers/python/index.html) (a C library) use **nanobind**
  and **scikit-build-core**, with stable-ABI wheels and zero-copy NumPy arrays.
- [Cantera](https://cantera.org/3.1/develop/distribution-packages/pypi-sdist-wheel.html) uses Cython, scikit-build-core and cibuildwheel;
  [RDKit](https://www.rdkit.org/docs/Install.html) is moving to nanobind; [python-igraph](https://github.com/igraph/python-igraph) uses a hand-written C-API
  extension and cibuildwheel.
- The [Scientific Python packaging guide](https://learn.scientific-python.org/development/guides/packaging-compiled/) recommends scikit-build-core for
  CMake projects.

The choices protolife makes from this are in [architecture.md](architecture.md).
