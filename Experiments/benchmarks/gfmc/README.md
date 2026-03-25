# GFMC Benchmarks

This benchmark suite turns the notebook-driven GFMC examples into a scriptable,
public-facing regression surface.

## Goals

- Cover the main GFMC capabilities present in `Experiments/systems/`
- Separate quick smoke tests from slower accuracy benchmarks
- Use theory-backed references whenever the system provides one
- Save CSV tables, figures, and debug logs under `outputs/`

## Tiers

- `smoke`: short runs with debug logging enabled by default. These are meant to
  catch broken code paths, invalid node handling, unstable branching, or bad
  warm-start plumbing quickly.
- `sweep`: moderate runs over a small time-step ladder. These are used to find
  stable plateaus and confirm expected trends before locking a final setup.
- `final`: the recommended benchmark-grade setup for each system. These runs
  emphasize lower time-step bias, healthier effective population, and figure
  generation.

## Covered Systems

- `free_particle_ring`: periodic wrapping and unguided generic GFMC
- `harmonic_oscillator_unguided`: open-boundary unguided propagation
- `harmonic_oscillator_guided`: exact importance sampling and warm start
- `cosine_lattice_ring`: guided vs unguided periodic external-potential test
- `hydrogen_unguided`: singular-potential stress test
- `hydrogen_fixed_node`: singular fixed-node stress test
- `two_particle_ho_guided`: exact two-coordinate guided benchmark
- `two_particle_ho_fixed_node`: exact odd-sector fixed-node benchmark
- `fermion_ring_fixed_node`: specialized fermion-ring kernel diagnostic benchmark
- `hardcore_boson_ring`: periodic pair-interaction stress benchmark
- `periodic_ion_single_particle`: exact-orbital warm-start benchmark
- `periodic_ion_spinless_fermions`: exact-determinant fixed-node benchmark
- `periodic_ion_bosons_jastrow`: bosonic orbital-Jastrow benchmark with energy
  brackets
- `periodic_ion_bosons_tg`: bosonic `|det|` Tonks-Girardeau scaffold benchmark

## Reference Policy

The suite uses three reference classes:

- Exact energy: free particle, harmonic oscillator, coupled HO, fixed-node odd
  HO, cosine lattice one-body problems, and determinant-based periodic-ion
  fermion benchmarks.
- Theory context: the bosonic periodic-ion cases report the noninteracting
  boson energy `N * ε0` and the Tonks-Girardeau scaffold energy `sum εα` as
  useful reference scales, but not as strict bounds because the smooth
  boson-boson repulsion can raise the exact energy above the pure TG scaffold.
- Stress-only diagnostics: the singular hydrogen and hard-core boson cases are
  tracked mainly through stability metrics, effective population, acceptance,
  and pair-separation diagnostics. The specialized fermion-ring kernel is
  currently also kept in this diagnostic bucket until its exact-reference
  behavior is tighter.

## Running

From the repository root:

```bash
julia --project=. Experiments/benchmarks/gfmc/run_gfmc_benchmarks.jl smoke all
julia --project=. Experiments/benchmarks/gfmc/run_gfmc_benchmarks.jl sweep accuracy
julia --project=. Experiments/benchmarks/gfmc/run_gfmc_benchmarks.jl final periodic_ion_single_particle,periodic_ion_spinless_fermions
julia --project=. Experiments/benchmarks/gfmc/generate_literature_report.jl
julia --project=. Experiments/benchmarks/gfmc/generate_theory_comparison_gallery.jl
```

Selections may be:

- `all`
- `accuracy`
- `stress`
- a comma-separated list of case ids

Outputs are written under `Experiments/benchmarks/gfmc/outputs/`:

- `tables/`: per-run histories and summary CSV files
- `figures/`: history, sweep, and diagnostic figures
- `logs/`: step-level debug logs for smoke runs or manually enabled debug runs

The literature-comparison report is written under
`Experiments/benchmarks/gfmc/outputs/literature/`:

- `tables/final_comparison.csv`: final paper-backed comparisons
- `tables/series_comparison.csv`: smoke, sweep, and final trajectories against
  the same theoretical references
- `figures/final_vs_literature.png`: final-value dumbbell plot
- `figures/final_error_vs_literature.png`: final-value error bars against paper
  references
- `figures/stress_paths_vs_literature.png`: time-step/error stress plot for the
  paper-backed cases
- `report.md`: markdown summary with theory notes and citation links

The source mapping from benchmark cases to primary references lives in
`Experiments/benchmarks/gfmc/LITERATURE_REFERENCES.md`.

For a tracked, GitHub-visible figure set instead of ignored benchmark outputs,
use `generate_theory_comparison_gallery.jl`. It reruns the paper-backed
benchmark cases, generates curated history/error overlays and density
comparisons, and writes them under `docs/assets/gfmc_theory_comparisons/` with
an index page at `docs/GFMC_BENCHMARK_THEORY_COMPARISONS.md`.

## Parameter Strategy

The tiered parameter choices follow a common pattern:

- Smoke tests use small walker populations and short imaginary time, but keep
  the time step conservative enough to expose drift, node, and branching bugs.
- Sweeps vary `dt` by factors of about two while keeping the total projection
  time roughly fixed for each case.
- Final runs choose the smallest stable `dt` in the sweep plateau and increase
  walker count to reduce noise.

The more singular or correlated the trial state is, the shorter the
reconfiguration interval and the smaller the recommended `dt`:

- Exact or near-exact guided cases can use `reconfiguration_interval = 5` or
  `10` with weaker feedback.
- Smooth unguided cases use moderate feedback and `reconfiguration_interval = 5`.
- Singular or hard-core cases tighten to `reconfiguration_interval = 1` or `2`
  and smaller `dt`.

## Notes

- The benchmark harness is intentionally script-first even when the original
  experiment is notebook-first. The notebooks remain the interactive analysis
  layer; this suite is the reproducible regression layer.
- Smoke runs enable the GFMC debug stream by default and write it to per-run
  log files instead of flooding stdout.
