# Experiments

This directory contains the curated, public-facing experiment layer for the
repository. Experiments are organized by physical system and then by Monte
Carlo method.

## Layout

```text
Experiments/
  benchmarks/              reproducible benchmark harnesses
  common/                  shared notebook helpers
  env/                     optional notebook environment
  systems/
    <system>/
      dmc/
        scripts/           stable CLI-style entrypoints, when available
        notebooks/         source notebooks for interactive analysis
      gfmc/
        notebooks/
      vmc/
        scripts/
        notebooks/
```

Not every system exposes every method, and not every method has a script
entrypoint. Notebook-only experiments are explicitly documented in
[`systems/README.md`](./systems/README.md).

## Policy

- Scripts are the preferred reproducible entrypoints when a scripted version exists.
- Notebooks are committed without saved outputs.
- Generated plots, CSV tables, animations, and frame dumps belong under
  `Experiments/**/outputs/`, which is gitignored.
- Figures should not be committed under source notebook directories.

## Running Scripts

Examples from the repository root:

```bash
julia --project=. Experiments/systems/free_particle_ring_1d/dmc/scripts/run_free_particle_ring_pbc.jl
julia --project=. Experiments/systems/cosine_lattice_ring_1d/dmc/scripts/run_cosine_lattice_ring_pbc.jl
VMC_TARGETN=2000 VMC_NSTEPS=40 VMC_WRITE_GIF=false \
  julia --project=. Experiments/systems/harmonic_oscillator_1d/vmc/scripts/vmc_animation_proposal_compare.jl
```

Many scripts accept environment-variable overrides. Check the configuration
block near the top of the script before launching large jobs.

The GFMC benchmark harness lives under [`benchmarks/gfmc/`](./benchmarks/gfmc/):

```bash
julia --project=. Experiments/benchmarks/gfmc/run_gfmc_benchmarks.jl smoke all
julia --project=. Experiments/benchmarks/gfmc/run_gfmc_benchmarks.jl sweep accuracy
julia --project=. Experiments/benchmarks/gfmc/run_gfmc_benchmarks.jl final stress
```

## Running Notebooks

Notebook setup cells search upward for the repository root and activate the
appropriate Julia project. If you want the optional notebook environment, run:

```bash
julia --project=Experiments/env -e 'using Pkg; Pkg.instantiate()'
```

Then open the notebook of interest in your preferred Julia notebook frontend.

## Experiment Index

The curated system inventory, along with the intended entrypoint for each
system, lives in [`Experiments/systems/README.md`](./systems/README.md).
