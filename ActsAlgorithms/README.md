# ActsAlgorithms — ALICE 3 add-on ACTS algorithms

This subproject contains ALICE-3-specific algorithms, generators and validation
tools built on top of the ACTS Examples framework, plus the Python bindings that
expose them. Everything lives in the `AliceActsTrk` C++ namespace and is reached
from Python through the `AliceActsPythonBindings` module. None of it modifies
ACTS itself — it is additive and links against the installed ACTS Examples
libraries.

Layout:

- **`ExampleAlgorithms/`** — the algorithms and generators (a shared library
  `ExampleAlgorithms`).
- **`Validation/`** — track-fit performance / residual-pull writers and tools.
- **`python/`** — the pybind11 module (`AliceActsPythonBindings`) that binds the
  above for use from ACTS Python job scripts.

---


### ExampleAlgorithm

A minimal skeleton `IAlgorithm` (single `testString` config). Template/starting
point for new algorithms; does no physics.

### TrackTruthMatcher

Hit-based track↔particle truth matching: for each track it counts the hits from
each truth particle (via `inputMeasurementParticlesMap`) and assigns the majority
particle above `matchingRatio`. Produces `track_particle_matching` /
`particle_track_matching` (with looper-protection options). Purely hit-based — no
time or position.

### HitRemoverAlgorithm

Rebuilds a measurement collection with a chosen subset of hits removed (e.g. hits
already used by accepted tracks), remapping the measurement/particle index maps
consistently. Used to feed later seeding/finding iterations with the remaining
hits.

### TrackMergerAlgorithm

Concatenates several input track collections (e.g. from multiple seeding/finding
iterations) into one output collection. Deduplication/clone removal is done by
the ambiguity-resolution step, not here.

### LutMaker

Builds look-up tables (e.g. per-`(eta, pT)` 2D histograms) from truth-matched
tracks for parameterised performance/resolution studies.

### BunchCrossingVertexGenerator

Vertex-position generator with a **discrete bunch-crossing time** (header
`AliceVertexGenerators.hpp`): `x, y, z` are Gaussian while the time is drawn
uniformly among `nBunchCrossings` slots spaced by `bunchSpacing`, centred on 0.
A drop-in replacement for `GaussianVertexGenerator` when simulating a bunch
train.

### PerPixelNoiseAlgorithm

Injects fake single-pixel noise so the reconstruction sees realistic detector
noise (fake hits/clusters → fake, contaminated tracks) that standard digitization
never produces — needed to study fake-hit rates, cluster sizes and tracking
purity.

Per sensitive pixel surface in a configurable `(volume, layer)` set it draws
`Poisson(nPixels * noiseProbPerPixel)` noise cells, then **re-clusters**: real
cells (from the input `ClusterContainer`) and noise cells go through the ACTS
`ModuleClusters` engine (the digitizer's own clustering), so a noise cell next to
a real cluster is absorbed and adjacent noise cells merge. Rebuilt clusters
inherit the truth of their real cells; **pure-noise clusters get no truth entry**
(genuine fakes, identifiable by an empty contributor list).

Runs after digitization and before space-point making. Rename the digitizer
outputs (e.g. `_raw`) and let it produce the canonical names, so nothing
downstream changes. Python wrapper: `addPerPixelNoise` (`alice3_perpixelnoise.py`).

Key config: `trackingGeometry`, `digitizationConfigs`, `randomNumbers`,
`noiseProbPerPixel` (default `1e-7`), `noiseActivation`, `layersByVolume`
(`{volume → [layers]}`; empty list = all layers, empty map = no-op),
`mergeCommonCorner`/`mergeNsigma`, and the `input*`/`output*` names.

### AliceVertexTruthMatcher

A copy of the ACTS `VertexTruthMatcher` (same logic and outputs, drop-in
replacement) with detailed per-event printouts and optional time-/primary-aware
matching — needed because the stock matcher is time-blind and, in overlapping
pile-up, mis-matches vertices and reports inflated contamination.

Printouts (per reco vertex): `x,y,z,t` ± uncertainties; the matched truth vertex
and Clean/Merged/Split class; and the track fractions that (1) match a truth
particle, (2) match the majority truth vertex, (3) come from other vertices
("contamination"), (4) fail the time-window cut. Each track shows a `fit:` line
(at the vertex — ACTS zeroes its per-track impact parameters and time) and an
`orig:` line (its own beamline perigee, with the real `z0`/`t0`). Truth vertices
are printed with their 4-position and contributing particles.

Optional modes (both default OFF → identical to `VertexTruthMatcher`):
`useTimeWindow`/`timeWindow` (count a track toward a truth vertex only if its
time is within `±timeWindow`; needs `inputTruthVertices`) and
`collapseSecondaryVertices` (key by the **primary** pp collision, so decay
daughters count toward their primary — for primary-vertex purity).


---

## Build / registration notes

- Sources are listed in `ExampleAlgorithms/CMakeLists.txt`. The library links the
  ACTS Examples `.so`s (`ActsExamplesFramework`, `ActsExamplesDigitization`) found
  via `CMAKE_PREFIX_PATH`. ACTS installs the Framework/EventData headers (seen via
  `Acts::Core`'s include dir) but **not** the Digitization headers used by
  `PerPixelNoiseAlgorithm`, so the CMake adds the Digitization *source* include
  dir (pass `-DACTS_SOURCE_DIR=/path/to/acts` if it is not auto-found).
- Python bindings are declared in `python/src/ExampleAlgorithmBindings.hpp` and
  registered in `python/src/AliceActsAlgorithmsModule.cpp`. Add a new component by
  writing its `addXxx` there and calling it in the module.
