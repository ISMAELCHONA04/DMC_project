# Experiment Systems

This index describes the systems kept in the public repository and the
recommended way to run each one.

For cross-system GFMC regression and benchmarking, prefer the scriptable suite
under `Experiments/benchmarks/gfmc/` instead of launching each notebook by
hand.

## Free Particle Ring 1D

- Purpose: minimal periodic-boundary smoke test for DMC and generic GFMC.
- Stable script: `free_particle_ring_1d/dmc/scripts/run_free_particle_ring_pbc.jl`
- Notebook entrypoints: DMC and GFMC notebooks under `free_particle_ring_1d/*/notebooks/`

## Cosine Lattice Ring 1D

- Purpose: one-particle periodic external-potential benchmark with guided vs unguided comparisons.
- Stable script: `cosine_lattice_ring_1d/dmc/scripts/run_cosine_lattice_ring_pbc.jl`
- Notebook entrypoints: DMC and GFMC notebooks under `cosine_lattice_ring_1d/*/notebooks/`

## Harmonic Oscillator 1D

- Purpose: simplest open-boundary benchmark across DMC, VMC, and GFMC.
- Stable scripts:
  - `harmonic_oscillator_1d/dmc/scripts/dmc_animation_single_particle.jl`
  - `harmonic_oscillator_1d/dmc/scripts/dmc_animation_single_particle_importance.jl`
  - `harmonic_oscillator_1d/vmc/scripts/vmc_animation_proposal_compare.jl`
- Notebook entrypoints: method notebooks under `harmonic_oscillator_1d/*/notebooks/`

## Hydrogen 1D

- Purpose: singular and regularized one-dimensional Coulomb-like trial-state examples.
- Stable scripts:
  - `hydrogen_1d/vmc/scripts/vmc_animation_hydrogen_proposal_compare.jl`
  - `hydrogen_1d/vmc/scripts/vmc_animation_regularized_h1d_proposal_compare.jl`
- Notebook entrypoints: DMC, VMC, and GFMC notebooks under `hydrogen_1d/*/notebooks/`
- Notes: the DMC and GFMC variants are notebook-driven in the current public tree.

## Two-Particle Harmonic Oscillator 1D

- Purpose: interacting two-particle comparison problems, including a fixed-node fermion case.
- Stable script entrypoints: notebook-driven in the current public tree.
- Notebook entrypoints: DMC and GFMC notebooks under `two_particle_ho_1d/*/notebooks/`

## Fermion Ring 1D

- Purpose: specialized spinless-fermion periodic-ring model with fixed-node DMC and GFMC.
- Stable scripts:
  - `fermion_ring_1d/dmc/scripts/run_spinless_fermion_ring_1d_pbc_fixed_node.jl`
  - `fermion_ring_1d/dmc/scripts/run_spinless_fermion_ring_1d_pbc_fixed_node_replicas.jl`
- Notebook entrypoints: DMC and GFMC notebooks under `fermion_ring_1d/*/notebooks/`
- Notes: the old TPU/JAX side path was removed from the public repository to keep the experiment surface maintainable.

## Hardcore Boson Ring 1D

- Purpose: bosonic Jastrow-guided GFMC notebook on a periodic cosine lattice with pair repulsion.
- Stable script entrypoints: notebook-driven in the current public tree.
- Notebook entrypoint: `hardcore_boson_ring_1d/gfmc/notebooks/HardcoreBosonRing1D_Coulomb_GFMC.ipynb`

## Periodic Ion Ring 1D

- Purpose: exact-orbital and bosonic/fermionic GFMC benchmarks on a ring with fixed ions.
- Stable script entrypoints: notebook-driven in the current public tree.
- Notebook entrypoints:
  - single-particle exact-orbital warm-start GFMC
  - spinless-fermion determinant warm-start GFMC
  - bosonic orbital-Jastrow warm-start GFMC
  - bosonic `|det|` Tonks-Girardeau scaffold warm-start GFMC
- Notes: this system is the most general ring benchmark in the repository and supports arbitrary ion positions through the helper layer in `periodic_ion_ring_helpers.jl`.
