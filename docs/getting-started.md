# Getting started

This guide takes you from a fresh machine to your own experiment: build protolife, run
a protocell next to its control, watch both, change the conditions, then try the larger
experiment with Genesis's sugar chemistry and use the library in your own code. Every
command below runs in the project's continuous integration on Linux, macOS and Windows.

## 1. What you need

- **A C compiler and CMake 3.24 or later.**
  - macOS: `xcode-select --install`, then `brew install cmake` (or `pip install cmake`).
  - Ubuntu or Debian: `sudo apt install build-essential git cmake`. If `cmake --version`
    is older than 3.24 (Ubuntu 22.04), use `pip install cmake`.
  - Windows: Visual Studio 2022 with "Desktop development with C++", which brings CMake;
    run the commands below in its "Developer PowerShell".
- **An internet connection for the first build**: SUNDIALS 7.9.0 (26 MB) is downloaded
  and built with the library.
- **Godot 4.7 or later** to watch the simulations: download it from
  [godotengine.org](https://godotengine.org/download). Nothing else is installed.

## 2. Build and test

```sh
git clone https://github.com/nicolashedoire/protolife.git
cd protolife
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure --parallel 4
```

The tests take about a minute. They include the reproduction of a published model
(Taneja and Higgs 2025), checks of every balance of matter, and a shortened step 8 of
Genesis: in one cycle a lineage whose chemistry runs takes a larger share of the
membranes, its twins without chemistry stay exactly even, and a lineage whose chemistry
stops gains nothing more. These three take most of that minute; `--label-exclude long`
leaves them out.

The programs are in `build/examples/` (on Windows `build\examples\Release\`, with
`.exe`). Run them from the project's folder: they write their recordings in
`recordings/`.

## 3. A first experiment: a protocell and its control

```sh
build/examples/first_experiment
```

It runs in a few seconds. Two protocells of Taneja and Higgs face the same conditions:
the active one holds an autocatalytic reaction (food makes catalyst, catalyst makes
more catalyst and decays into waste); its control holds the same food and no catalyst.
The table it prints shows the active cell swelling, stretching and dividing again and
again, while the control's volume stays put and its membrane keeps growing until it
becomes a thin tube, where the model ends, as its authors' program does.

Watch both:

```sh
godot --path viewer
```

(On macOS, if `godot` is not on your path: `/Applications/Godot.app/Contents/MacOS/Godot --path viewer`.)
The viewer lists the recordings; open `first-active`, then `first-control`. The membrane
has the shape the model computes, a sphere or a capsule; the points inside are the
catalyst (amber), the food (blue) and the waste (grey); at a division the daughter let
go drifts away. Play, pause, scrub and change the speed with the bar at the bottom; drag
to turn around the cell.

## 4. Change the conditions

Open `examples/first_experiment.c`. Its settings are at the top:

| Setting | What it is | Try |
|---|---|---|
| `LIPID_OUTSIDE` | lipid in the solution, in units of the balance C* | `1.0`: no growth of the membrane, the cell stops at 27/8 of its volume |
| `FLIP_FLOP` | lipid flipping between the membrane's two leaflets | `0.2`: no division at all; `0.01`: very unequal daughters |
| `WASTE_LEAVES` | how fast the waste leaves the cell | `0.01`: no steady size, the cell swells without end |
| `DURATION` | time simulated | longer runs show more cycles |

Rebuild and run again:

```sh
cmake --build build --config Release && build/examples/first_experiment
```

Each of these outcomes is a result of Taneja and Higgs's paper, and the test suite checks
them (`tests/test_capsule.c`).

## 5. The larger experiment: sugar chemistry inside dividing vesicles

`reproduction` puts an organisation of formose chemistry (formaldehyde becoming sugars,
from Genesis's Origin) inside 500 nm oleate vesicles in a dialysis bag, then runs cycles
of dilution, feeding with micelles and shear. With `living=0` no vesicle reacts: the
control.

```sh
build/examples/reproduction out=recordings/step-6-random2                 # about 8 minutes
build/examples/reproduction out=recordings/step-6-control living=0        # about 1 minute
```

In the living run each seeded vesicle becomes 16 in four cycles and every daughter
remakes its sugars; in the control the sugars are only shared out, down to about a
thousandth. At the end each run prints its balances of matter (fatty acid, buffer,
carbon) and how far they drifted: below a millionth of a millionth.

Its options are written `key=value`:

| Option | Default | Meaning |
|---|---|---|
| `network=` | `data/networks/e0-random2.net` | the chemistry inside: the five E0 cores are in `data/networks/` |
| `living=` | `1` | `0` for the control |
| `radius=` | `500` | radius of the vesicles, nm |
| `seeded=` | `0.9` | share of the fatty acid in vesicles holding the chemistry |
| `food=` | `0.5` | formaldehyde in the fresh medium, mol/L |
| `cycles=`, `period=` | `4`, `900` | cycles of dilution, feeding and shear, and their length (s) |
| `micelles=`, `feeds=` | `1`, `1` | equivalents of membrane fed per cycle, and in how many doses |
| `stirred=`, `shear=` | `1`, `20` | shear all along (1) or at two thirds of each cycle (0), and its rate (1/s) |
| `start=`, `dilution=` | `3600`, `0.5` | installation before the first cycle (s), and share of the bag renewed per cycle |
| `seed=` | `core` | `fresh`: the seeded vesicles hold fresh medium too, the twins of the others (a neutral control) |
| `stop=` | `0` | `1`: the seeded vesicles' chemistry stops at the first cycle; they keep what it made |
| `vesicles=`, `draw=` | `0`, `1` | a finite bag holding that many seeded vesicles, its dilutions drawn at random with that seed; `0` follows the mean |

## 6. Your own chemistry

A network is a text file; this is `data/networks/taneja-higgs.net`:

```
species C1 C 0 dynamic 0.5
species C2 C2 0 dynamic 0.5
species W C 0 dynamic 0
reaction 2 C1 + 2 C2 -> 3 C2 ; k=1 ; source=Taneja and Higgs 2025, k
reaction C2 -> 2 W ; k=0.01 ; source=Taneja and Higgs 2025, w
```

A species has a name, a formula, a charge, a role (`dynamic`, `reservoir` held fixed
from outside, or `implicit` folded into the rate constants) and a concentration in
mol/L. A reaction has an equation, a mass-action rate constant and a source, which is
required. protolife checks that every reaction keeps its atoms and its charge, and
refuses the file otherwise, naming the reaction. The full format is described in
`include/protolife/network.h`.

## 7. Use protolife in your own program

Install it as a shared library, which carries everything it needs:

```sh
cmake -S . -B build-shared -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
cmake --build build-shared --config Release
cmake --install build-shared --config Release --prefix $HOME/.local
```

Then, in your project's `CMakeLists.txt`:

```cmake
find_package(protolife REQUIRED)    # pass -DCMAKE_PREFIX_PATH=$HOME/.local if needed
add_executable(my_experiment main.c)
target_link_libraries(my_experiment PRIVATE protolife::protolife)
target_compile_definitions(my_experiment PRIVATE NETWORKS="${protolife_DATA_DIR}/networks")
```

`tests/install/` is such a project; it runs in the continuous integration.

## 8. Reading the results

- Every call that can fail returns a status, and its diagnostic says what went wrong,
  at what time, and in which quantity. A failed call leaves the simulation as it was.
- Simulations of vesicles and reactors check their balances of matter at every step and
  refuse a step that breaks one.
- Each mechanism names its source in the headers (`include/protolife/`); the parameters
  that are assumed rather than measured are said to be.

Questions and reports are welcome in the
[project's issues](https://github.com/nicolashedoire/protolife/issues).
