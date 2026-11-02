# protolife

![Illustration of protocell vesicles and prebiotic chemistry](docs/images/protolife-cover.png)

**Building blocks for synthetic life, in simulation.** protolife is meant to be to
protocells what a neural-network library is to deep learning: a set of composable,
tested parts (chemistry, compartments, media, protocols, measurements) that anyone
can assemble to build a candidate living system and test it against a control.

> **Status: pre-alpha.** The library is being extracted from Genesis, a C11 research
> code that already simulates self-maintaining sugar chemistry inside fatty-acid
> vesicles, their growth, division and one-bit heredity. The C core has started with
> its foundation: error reporting and stiff integration over SUNDIALS. The design is in
> [the landscape of existing tools](docs/landscape.md) and
> [the architecture](docs/architecture.md). The Python API shown below does not exist yet.

## What it is for

Laboratories can make the building blocks of life, compartments that grow and split,
and chemistries that amplify themselves. Putting them together into something that
persists *and passes its organisation on* is the open problem. protolife aims to let
anyone explore that problem in physical units, with the bookkeeping that keeps a
simulation honest.

No existing open library combines these in one engine (see
[docs/landscape.md](docs/landscape.md)):

- vesicle mechanics in SI units: osmosis, membrane tension, rupture, growth by
  micelle feeding, shear-driven division;
- a population of compartments with lineage tracking;
- a stiff chemistry running inside each compartment, deterministic or stochastic;
- chemistry generated from molecular graphs and reaction rules, with
  thermochemistry and kinetics;
- laboratory protocols: feeding, dilution and shear cycles, dialysis;
- mechanisms calibrated on published experiments, each with its validity domain.

## Principles

1. **Every mechanism cites its source** and declares the range where it holds.
   Parameters outside that range are refused, not extrapolated in silence.
2. **Conservation is checked.** Matter balances close at every step, and the drift is
   reported.
3. **A solver that fails says so.** It never returns a silently corrected state.
4. **Controls are built in.** Every experiment can be rerun without the mechanism
   under test (a lifeless twin, a random-parent baseline).
5. **Nothing is programmed to succeed.** No cell, division, heredity or adaptation is
   written as an outcome. A model may produce nothing, and that is a result.

## Planned API

A sketch of the Python layer, not implemented yet: two lineages of vesicles in one
dialysis bag, one carrying a self-maintaining chemistry and one without it, through
eight cycles of feeding and shear.

```python
import protolife as pl

bag = pl.Medium.dialysis_bag(food=pl.formaldehyde(0.5))       # mol/L
living = pl.Vesicles.oleate(radius=600e-9, chemistry=pl.formose.core("carbons"))
unlit = pl.Vesicles.oleate(radius=420e-9, chemistry=None)

world = pl.World(bag, [living, unlit])
trace = world.run(pl.Cycles(count=8, period=15 * 60, shear="permanent"))
trace.plot("sugar_carbon", by="lineage")
```

## Building the C core

Requirements: CMake 3.24 or later and a C11 compiler. SUNDIALS 7.9.0 is downloaded and
built during configuration (to build offline, pass
`-DFETCHCONTENT_SOURCE_DIR_SUNDIALS=/path/to/sundials-7.9.0`).

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

What exists today, in C:

- `protolife/diagnostic.h`: every call that can fail returns a status and fills a
  diagnostic with the code of the library underneath, the model time and a message.
- `protolife/ode.h`: stiff integration by CVODE with event functions. A call that cannot
  reach its time leaves the state as it was and says why; nothing is clipped. On
  Robertson's stiff chemistry it matches the Genesis integrator it replaces to within
  2.4e-9, in about 80 times fewer steps.

## Roadmap

1. **Landscape of existing tools**: done, see [docs/landscape.md](docs/landscape.md).
2. **C11 core** (in progress): the general bricks behind a stable C API, with
   regression tests against the Genesis results they come from. Done: diagnostics,
   stiff integration. Next: reaction networks, stochastic simulation, molecules and
   rules, vesicles, media and protocols, observation.
3. **Python layer**: `pip install protolife`, with the example above as the first
   working script.
4. **Documentation and a catalogue of calibrated worlds**, plus reference tests
   against published protocell models.

## License

BSD 3-Clause, see [LICENSE](LICENSE).
