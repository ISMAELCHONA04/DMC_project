# VMC User Guide

## Overview

The VMC implementation lives under `src/UseCases/vmc/` and supports:

- Metropolis sampling with a trial wave function
- drifted-Gaussian proposals through `DriftGaussianProposal()`
- plain Gaussian proposals through `GaussianProposal()`
- open or periodic boundaries inherited from the Hamiltonian
- raw position inputs or explicit `Walker` objects
- exporting final equilibrated positions through `sample_vmc_positions(...)` for warm-starting GFMC

Unlike DMC, the current VMC workflow does not use a separate node policy.
Nodal structure only enters indirectly through the trial wave function if you
choose to encode sign information for analysis.

## Main Types

### `VMCParams`

Use either the positional constructor

```julia
VMCParams(dt, nsteps, targetN, ET0)
```

or the keyword constructor

```julia
VMCParams(; dt, nsteps, targetN, ET0=0.0)
```

Parameters:

- `dt`: proposal time step
- `nsteps`: number of Metropolis sweeps
- `targetN`: number of walkers/configurations to evolve
- `ET0`: optional reference value retained for compatibility with existing scripts and plot annotations

Validation rules:

- `dt > 0`
- `nsteps >= 1`
- `targetN >= 1`

Important note:

- the current core VMC algorithm does not dynamically use `ET0`
- experiment scripts sometimes keep it as a convenient reference energy for labeling or comparisons

### `TrialWF`

VMC uses the same trial-wavefunction interface as DMC:

```julia
TrialWF(logpsi, gradlogpsi, lapllogpsi)
TrialWF(logpsi, gradlogpsi, lapllogpsi, signpsi)
```

Callbacks:

- `logpsi(R)`: returns `log(abs(psi_T(R)))`
- `gradlogpsi(R)`: returns the gradient of `log(abs(psi_T))`
- `lapllogpsi(R)`: returns the Laplacian of `log(abs(psi_T))`
- `signpsi(R)`: optional sign callback

The VMC code uses:

- `logpsi` in the Metropolis acceptance ratio
- `gradlogpsi` in drifted-Gaussian proposals
- `lapllogpsi` and `gradlogpsi` in the local-energy estimator

### Proposal modes

Available proposal policies:

- `DriftGaussianProposal()`
- `GaussianProposal()`

Symbol aliases:

- `:drift_gaussian`, `:drifted_gaussian`, `:drift`
- `:gaussian`, `:position_gaussian`

`DriftGaussianProposal()` is usually the better default when the trial state is
smooth and informative. `GaussianProposal()` is useful as a simpler baseline.

## Defining the Physics Input

### Hamiltonian and boundaries

VMC uses the same Hamiltonian and boundary-condition objects as DMC:

```julia
Hamiltonian(N, D, V, bc=OpenBoundary1D())
```

with:

- `N`: number of coordinates in the configuration vector
- `D`: diffusion constant used by the proposal kernel
- `V`: potential callback `R -> scalar`
- `bc`: boundary-condition policy

Examples:

```julia
H_open = Hamiltonian(1, 0.5, R -> 0.5 * R[1]^2)
H_ring = Hamiltonian(1, 0.5, R -> cos(2pi * R[1]), PeriodicBoundary1D(0.0, 1.0))
```

For periodic systems, initial configurations and proposed moves are wrapped
back into the simulation cell automatically.

### Trial-wavefunction examples

Gaussian harmonic-oscillator ansatz:

```julia
trial = TrialWF(
    R -> -0.5 * R[1]^2,
    R -> [-R[1]],
    R -> -1.0,
)
```

Two-coordinate separable ansatz:

```julia
trial = TrialWF(
    R -> -0.25 * (R[1]^2 + R[2]^2),
    R -> [-0.5 * R[1], -0.5 * R[2]],
    R -> -1.0,
)
```

If you want to track sign information for an odd state, you can still provide a
fourth callback:

```julia
trial = TrialWF(
    R -> log(max(abs(R[1]), 1e-8)) - 0.5 * R[1]^2,
    R -> abs(R[1]) < 1e-8 ? [0.0] : [inv(R[1]) - R[1]],
    R -> abs(R[1]) < 1e-8 ? -1.0 : (-(inv(R[1])^2) - 1.0),
    R -> abs(R[1]) < 1e-8 ? 0.0 : sign(R[1]),
)
```

The current VMC implementation does not reject node-crossing proposals, but the
same `TrialWF` object can be reused later in DMC with `FixedNode()`.

## Running a Simulation

### Recommended one-call interface

```julia
using Random
using System1D

H = Hamiltonian(1, 0.5, R -> 0.5 * R[1]^2)
trial = TrialWF(
    R -> -0.5 * R[1]^2,
    R -> [-R[1]],
    R -> -1.0,
)
params = VMCParams(; dt=0.02, nsteps=200, targetN=4096, ET0=0.5)
positions = [[randn()] for _ in 1:params.targetN]

sim = run_vmc(H, params, positions, trial;
    rng=MersenneTwister(7),
    proposal=:drift_gaussian,
    snapshot_steps=[0, params.nsteps],
    show_progress=true,
)
```

### Sampling Warm-Start Positions For GFMC

If you only want a VMC-equilibrated ensemble to hand off into GFMC, use
`sample_vmc_positions(...)` instead of `run_vmc(...)`:

```julia
using Random
using System1D

H = Hamiltonian(1, 0.5, R -> 0.5 * R[1]^2, PeriodicBoundary1D(0.0, 1.0))
trial = TrialWF(
    R -> -0.5 * R[1]^2,
    R -> [-R[1]],
    R -> -1.0,
)
params = VMCParams(; dt=0.02, nsteps=100, targetN=256, ET0=0.0)
seed_positions = [[rand()] for _ in 1:params.targetN]

positions = sample_vmc_positions(
    H,
    params,
    seed_positions,
    trial;
    rng=MersenneTwister(11),
    proposal=:drift_gaussian,
)
```

This helper runs the VMC sweeps but only returns the final raw positions. It
does not keep the full walker-position history unless you explicitly run
`run_vmc(...)` yourself.

`run_vmc(...)` accepts:

- a vector of raw position vectors, such as `[[x1], [x2], ...]`
- a single raw position vector, for a one-walker simulation
- a matrix whose columns are walker configurations
- an explicit vector of `Walker` objects

Keyword arguments:

- `rng`
- `proposal`
- `step`
- `snapshot_steps`
- `show_progress`
- `progress_every`
- `progress_label`
- `debug`
- `debug_every`
- `debug_io`

### Low-level interface

If you want to construct the simulation object yourself:

```julia
sim = VMCSim(H, params, positions, trial, MersenneTwister(7); proposal=DriftGaussianProposal())
run_vmc!(sim; snapshot_steps=[0, 50, 100])
```

## Important Functions

### `VMCSim(...)`

Constructs the simulation state, wraps initial configurations into the
Hamiltonian cell, and builds an `ImportanceGuiding(trial, H)` object
internally.

### `initialize!(sim, configs)`

Resets an existing `VMCSim` with new walkers or raw configurations and clears:

- `energy_history`
- `energy_variance_history`
- `walker_positions_history`
- `acceptance_rate`
- `acceptance_count`

### `proposal_step!(sim, walker)`

Generates a proposal and its Metropolis log-ratio.

Behavior by proposal type:

- `DriftGaussianProposal()`: centered at `R + D * dt * drift`
- `GaussianProposal()`: centered at `R`

In both cases, periodic boundaries are handled automatically.

### `metropolis_step!(sim, walker)`

Applies the Metropolis accept/reject test and returns the accepted walker or
the original walker.

### `vmc_step!(sim)`

Runs one full sweep over the walker ensemble, updates `sim.step`, and refreshes
the cumulative `acceptance_rate`.

### `run_vmc!(sim; ...)`

Runs `sim.params.nsteps` sweeps and records:

- `energy_history`
- `energy_variance_history`
- `walker_positions_history`

Debug mode:

- set `debug=true` to print a per-step line with timing, walker count,
  acceptance information, and current energy statistics
- use `debug_every` to print every `n`th step instead of every step
- use `debug_io` if you want to capture the diagnostics instead of printing to stdout

### `sample_vmc_positions(H, params, configs, trial_wf; ...)`

This convenience helper is meant for cross-method handoff:

1. it builds a `VMCSim`
2. runs `run_vmc!` for `params.nsteps`
3. returns deep copies of the final walker positions as `Vector{Vector{Float64}}`

The accepted inputs and proposal keywords match `run_vmc(...)`.
Use this when you need a `|psi_T|^2`-distributed initial ensemble for GFMC but
do not need the full VMC history object.

### `estimate_energy_vmc(sim)` and `compute_variance(sim)`

Convenience wrappers around the shared energy estimators:

- `estimate_energy_vmc(sim)` returns the current mean local energy
- `compute_variance(sim)` returns the current local-energy variance

## Diagnostics and Outputs

The main quantities stored on a completed `VMCSim` are:

- `sim.energy_history`
- `sim.energy_variance_history`
- `sim.walker_positions_history`
- `sim.acceptance_rate`
- `sim.acceptance_count`

Typical workflow:

1. define `Hamiltonian` and `TrialWF`
2. build `VMCParams`
3. prepare initial configurations
4. run `run_vmc(...)`
5. inspect energy convergence, variance convergence, and accepted distributions

For projector warm starts, replace step 4 with `sample_vmc_positions(...)` and
feed the returned raw positions into `GFMCSim(...)` or
`run_gfmc_with_vmc_init(...)`.
