# System1D

`System1D` is a Julia repository for one-dimensional quantum Monte Carlo
experiments. It provides a small reusable package layer for defining
Hamiltonians, boundaries, trial states, and walkers, then builds DMC, VMC, and
fixed-population GFMC workflows on top of that package.

The repository is organized as a research codebase with two goals:

- a documented core API for 1D Monte Carlo methods
- a curated set of reproducible experiments under [`Experiments/`](./Experiments)

## Scope

The current public surface includes:

- diffusion Monte Carlo with optional importance sampling and fixed-node rejection
- variational Monte Carlo with drifted-Gaussian and symmetric Gaussian proposals
- fixed-population GFMC with reconfiguration and optional VMC warm starts
- open, periodic, and twisted-periodic one-dimensional boundaries
- a built-in `SpinlessFermionRing1D` model used by the fermion-ring workflows

GFMC in this repository is a weighted, fixed-population projector method for
continuous coordinates. It is documented and tested, but still narrower in
scope than the DMC and VMC layers.

## Repository Layout

```text
src/
  Domain/                  core physics types and geometry helpers
  UseCases/common/         shared utilities used by multiple methods
  UseCases/dmc/            diffusion Monte Carlo implementation
  UseCases/vmc/            variational Monte Carlo implementation
  UseCases/gfmc/           fixed-population GFMC implementation
test/
  runtests.jl              package-level regression tests
docs/
  *_USER_GUIDE.md          method guides
  GFMC_CODE_WALKTHROUGH.tex
Experiments/
  common/                  notebook helpers shared across experiments
  env/                     optional notebook environment
  systems/                 curated example systems and entrypoints
```

Generated outputs are intentionally excluded from version control. Scripts and
notebooks write figures, CSV tables, animations, and other byproducts into
`Experiments/**/outputs/`, which is gitignored.

## Installation

Instantiate the package environment from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the package tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

If you want to run the notebooks with optional extras such as notebook-side
plotting packages, instantiate the optional experiment environment as well:

```bash
julia --project=Experiments/env -e 'using Pkg; Pkg.instantiate()'
```

## Main Workflows

### Package usage

The package layer is intended to be usable directly:

```julia
using Random
using System1D

bc = PeriodicBoundary1D(0.0, 10.0)
H = Hamiltonian(1, 0.5, R -> 0.5 * R[1]^2, bc)
trial = TrialWF(
    R -> -0.5 * R[1]^2,
    R -> [-R[1]],
    R -> -1.0,
)

vmc_params = VMCParams(; dt=0.02, nsteps=100, targetN=256, ET0=0.5)
gfmc_params = GFMCParams(0.02, 100, 20, 256, 0.5, 0.2, 5, 5.0, 20)
seed_positions = [[rand()] for _ in 1:gfmc_params.targetN]

sim = run_gfmc_with_vmc_init(
    H,
    gfmc_params,
    seed_positions,
    trial,
    vmc_params;
    vmc_rng=MersenneTwister(11),
    gfmc_rng=MersenneTwister(12),
    proposal=:drift_gaussian,
    snapshot_steps=[0, gfmc_params.nsteps],
)
```

See the method guides in [`docs/`](./docs) for the full API surface.

### Reproducible experiments

Representative experiment entrypoints from the repository root:

```bash
julia --project=. Experiments/systems/free_particle_ring_1d/dmc/scripts/run_free_particle_ring_pbc.jl
julia --project=. Experiments/systems/cosine_lattice_ring_1d/dmc/scripts/run_cosine_lattice_ring_pbc.jl
VMC_TARGETN=2000 VMC_NSTEPS=40 VMC_WRITE_GIF=false \
  julia --project=. Experiments/systems/harmonic_oscillator_1d/vmc/scripts/vmc_animation_proposal_compare.jl
```

For the curated notebook and script inventory, see:

- [`Experiments/README.md`](./Experiments/README.md)
- [`Experiments/systems/README.md`](./Experiments/systems/README.md)

## Documentation Index

- [`docs/DMC_USER_GUIDE.md`](./docs/DMC_USER_GUIDE.md)
- [`docs/VMC_USER_GUIDE.md`](./docs/VMC_USER_GUIDE.md)
- [`docs/GFMC_USER_GUIDE.md`](./docs/GFMC_USER_GUIDE.md)
- [`docs/README.md`](./docs/README.md)
- [`CONTRIBUTING.md`](./CONTRIBUTING.md)

## Notes for Collaborators

- Notebooks are committed as source-only documents. Saved cell outputs should
  not be committed.
- Generated figures and tables belong under gitignored `outputs/` directories,
  not next to source notebooks.
- The repository intentionally excludes machine-specific Colab/TPU branches and
  other unmaintained side paths so the public tree stays coherent.
