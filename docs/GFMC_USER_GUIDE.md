# GFMC User Guide

## Overview

The GFMC code lives under `src/UseCases/gfmc/` and implements a
fixed-population projector Monte Carlo workflow.

In this repository, "GFMC" means:

- continuous-coordinate walkers, stored in a batched matrix
- Gaussian or drift-diffusion proposals in imaginary time
- multiplicative walker weights instead of explicit birth/death every step
- periodic reconfiguration back to exactly `targetN` walkers
- optional importance sampling and optional fixed-node handling (for fermions)

That makes this code closer to a weighted, fixed-population DMC variant than to
an explicit matrix-element lattice GFMC implementation. The important design
difference from the DMC code in this repository is:

- DMC changes the walker count through branching
- GFMC keeps the walker count fixed and instead tracks a weight per walker
- GFMC periodically resamples those weights into a new equally weighted ensemble

The main public entry points are:

- `GFMCParams`
- `GFMCSim`
- `run_gfmc!`
- `run_gfmc_with_vmc_init`

There is still no generic one-call `run_gfmc(...)` convenience wrapper, so the
normal low-level usage pattern is:

1. build a Hamiltonian or `SpinlessFermionRing1D` model
2. construct `GFMCParams`
3. construct `GFMCSim(...)`
4. call `run_gfmc!(sim; ...)`

For VMC-backed warm starts, there is now a dedicated convenience wrapper:

```julia
run_gfmc_with_vmc_init(...)
```

## Source Map

The implementation is split across these files:

- `src/UseCases/gfmc/abstract.jl`
  Defines the GFMC extension points: kernels, backends, and reconfiguration
  policies.
- `src/UseCases/gfmc/params.jl`
  Defines `GFMCParams`.
- `src/UseCases/gfmc/kernels.jl`
  Defines batched kernel evaluation, including the generic Hamiltonian-based
  kernel and the specialized `SpinlessFermionRing1DKernel`.
- `src/UseCases/gfmc/sim.jl`
  Defines `GFMCSim`, constructors, estimators, and state/history recording.
- `src/UseCases/gfmc/reconfiguration.jl`
  Defines the fixed-population resampling schemes.
- `src/UseCases/gfmc/run.jl`
  Defines the propagation step, reference-energy update, and main run loop.

The GFMC layer also depends on shared domain-level types:

- `src/Domain/hamiltonian.jl`
- `src/Domain/boundaries1D.jl`
- `src/Domain/trialwf.jl`
- `src/Domain/guiding.jl`
- `src/Domain/nodepolicy.jl`
- `src/Domain/spinless_fermion_ring_1d.jl`

## Main Types

### `GFMCParams`

Constructor:

```julia
GFMCParams(dt, nsteps, nequil, targetN, ET0, feedback, reconfiguration_interval, branch_cap, energy_window)
```

Fields:

- `dt`: imaginary-time step size
- `nsteps`: number of propagation steps
- `nequil`: number of equilibration steps before the code starts updating the
  reference energy `ET`
- `targetN`: fixed walker count
- `ET0`: initial reference-energy guess
- `feedback`: strength of the population-control feedback term used in
  `_update_reference_energy!`
- `reconfiguration_interval`: how many propagation steps occur between weight
  resampling events
- `branch_cap`: upper bound applied to a single-step branching weight
- `energy_window`: moving-average window used when smoothing the energy input
  to the `ET` update

Validation rules:

- `dt > 0`
- `nsteps >= 1`
- `0 <= nequil <= nsteps`
- `targetN >= 1`
- `reconfiguration_interval >= 1`
- `branch_cap > 0`
- `energy_window >= 1`

### `AbstractGFMCKernel`

The kernel is the physics adapter used by the GFMC engine. It tells the engine:

- how many coordinates each walker has
- what diffusion constant to use
- what boundary condition to use
- whether importance sampling is active
- whether fixed-node logic is active
- how to evaluate drift, `logpsi`, branch energy, and sign for a whole walker
  batch

Built-in implementations:

- `GenericGFMCKernel`
- `SpinlessFermionRing1DKernel`

### `GFMCBatchData`

`GFMCBatchData` is the per-walker auxiliary state stored alongside the position
matrix:

- `drift`: matrix of quantum-force values, same shape as the position matrix
- `logpsi`: vector of `log(abs(psi_T))`
- `branch_energy`: vector used in the branching-weight formula
- `sign`: vector used by fixed-node killing

There are always two copies:

- `state_data` for the current accepted configuration set
- `proposal_data` for the current proposed configuration set

### `GFMCSim`

`GFMCSim` is the full simulation state. The most important fields are:

- `positions`: current walker matrix with shape `(ncoords, nwalkers)`
- `proposal_positions`: proposal buffer with the same shape
- `state_data` and `proposal_data`: double-buffered `GFMCBatchData`
- `weights`: current walker weights
- `kernel`: physics adapter
- `params`: runtime parameters
- `rng`: random-number generator
- `backend`: storage and copy backend, currently `CPUBackend()`
- `reconfiguration`: `SystematicReconfiguration()` or
  `MultinomialReconfiguration()`
- `step`: current step counter
- `ET`: current reference energy

History arrays:

- `ET_history`
- `population_history`
- `energy_mean_history`
- `energy_variance_history`
- `effective_population_history`
- `mean_weight_history`
- `acceptance_history`
- `walker_positions_history`

A few implementation details matter:

- each column of `positions` is one walker configuration
- the walker count never changes, so `size(sim.positions, 2)` is always
  `targetN`
- the code uses buffer swapping instead of allocating a fresh proposal array
  every step

## Construction Flow

### Accepted position formats

The low-level constructor normalizes input into a `Matrix{Float64}` with one
walker per column. It accepts:

- a matrix whose columns are walker configurations
- a vector of position vectors
- a vector of `Walker` objects

When `Walker` objects are passed in, GFMC keeps only `position(w)`. Any stored
walker phase is discarded by the matrix conversion step.

Unlike the DMC/VMC helpers, there is no special single-vector convenience path
for a one-walker GFMC run. A single walker should therefore be passed as either
`[[x]]`, `reshape([x], 1, 1)`, or `[Walker([x])]`.

### Constructor sequence

`GFMCSim(kernel, params, positions, rng; ...)` does the following:

1. convert the input positions to a dense `Float64` matrix
2. verify that the number of columns equals `params.targetN`
3. verify that the number of rows matches `configuration_dimension(kernel)`
4. wrap all coordinates into the boundary cell with `wrap_configurations!`
5. allocate `proposal_positions`
6. allocate `state_data` and `proposal_data`
7. evaluate the current batch into `state_data`
8. allocate `weights` and set them all to `1.0`
9. initialize `ET = params.ET0`
10. clear all history vectors

There are also two higher-level constructor families:

```julia
GFMCSim(H, params, walkers_or_positions, rng;
    guiding=NoGuiding(),
    nodepolicy=NoNode(),
    backend=CPUBackend(),
    reconfiguration=SystematicReconfiguration(),
)
```

This path wraps the Hamiltonian into a `GenericGFMCKernel`.

```julia
GFMCSim(model::SpinlessFermionRing1D, params, positions, rng;
    use_guiding=true,
    nodepolicy=FixedNode(),
    backend=CPUBackend(),
    reconfiguration=SystematicReconfiguration(),
)
```

This path uses the specialized ring kernel.

### VMC-backed warm starts

If you want the initial GFMC walkers to be sampled from `|psi_T|^2` rather
than from a uniform or ad hoc seed distribution, use either:

- `sample_vmc_positions(H, vmc_params, init_configs, trial_wf; ...)`
- `run_gfmc_with_vmc_init(...)`

`sample_vmc_positions(...)` runs a short VMC equilibration and returns the
final walker positions as raw configuration vectors. Those positions can be
passed directly into `GFMCSim(...)`.

`run_gfmc_with_vmc_init(...)` is the one-call version that performs both the
VMC warm start and the GFMC projection.

## Kernel Evaluation

### Generic Hamiltonian path

`GenericGFMCKernel` is the general adapter for an arbitrary `Hamiltonian`.

With importance sampling enabled:

- `drift = 2 * gradlogpsi(R)`
- `logpsi = log(abs(psi_T(R)))`
- `branch_energy = local_energy(R)`
- `sign = signpsi(R)`

The local energy comes from `ImportanceGuiding`:

```julia
EL(R) = -D * (lapllogpsi(R) + sum(abs2, gradlogpsi(R))) + V(R)
```

With importance sampling disabled:

- `drift` is filled with zeros
- `logpsi` is set to `0.0`
- `branch_energy` is the bare potential `V(R)`
- `sign` is set to `1.0`

That last point is important: in the current implementation, the unguided
generic path does not carry nodal information.

### Specialized fermion-ring path

`SpinlessFermionRing1DKernel` is an optimized batch evaluator for
`SpinlessFermionRing1D`.

When `use_guiding=true`, the kernel computes in one pass:

- the one-body lattice term `lambda * cos(k_lat * x_i)` in `logpsi`
- the pairwise fermionic term `log(abs(sin(alpha_pair * dx_ij)))`
- the drift from the derivative of those terms
- the Laplacian contribution needed for the local energy
- the trial-wavefunction sign from the sign of the pairwise sine product

Its branch energy is the same local-energy formula used elsewhere:

```julia
EL(R) = -D * (lapl_logpsi(R) + |grad_logpsi(R)|^2) + V(R)
```

When `use_guiding=false`, it falls back to:

- zero drift
- zero `logpsi`
- branch energy equal to the cosine-lattice potential
- sign fixed to `+1`

As with the generic kernel, fixed-node logic only has useful sign information
when the guiding path is active.

## Propagation Step

The core propagation routine is `step!(sim)`.

Let:

- `D = diffusion_constant(sim.kernel)`
- `dt = sim.params.dt`
- `sigma = sqrt(2 * D * dt)`

### 1. Build proposals

The code first fills `sim.proposal_positions` with standard normal noise.

If importance sampling is enabled, each walker coordinate is updated as:

```julia
R' = R + D * dt * F(R) + sigma * eta
```

where `F(R)` is the drift stored in `sim.state_data.drift` and `eta` is a
standard normal variate.

If importance sampling is disabled, the proposal is plain diffusion:

```julia
R' = R + sigma * eta
```

After proposal generation, `wrap_configurations!` maps every coordinate back
into the boundary cell.

### 2. Evaluate proposed walkers

The proposed walker matrix is passed through `evaluate_configuration_data!`,
which fills:

- proposal drift
- proposal `logpsi`
- proposal branch energy
- proposal sign

These values are used both for acceptance and for weight updates.

### 3. Fixed-node killing

For each walker, the code checks fixed-node viability before any
Metropolis-style acceptance step:

- `sold = sim.state_data.sign[j]`
- `snew = sim.proposal_data.sign[j]`

The walker is killed if:

- `sold == 0.0`
- `snew == 0.0`
- `sold != snew`

When a walker is killed:

- the proposal position/data are overwritten with the old state
- `sim.weights[j] = 0.0`
- the loop continues to the next walker

This means fixed-node violations do not move the walker, but they do remove its
future influence until the next reconfiguration event.

### 4. Metropolis correction for guided proposals

If importance sampling is enabled, the code applies a Metropolis-Hastings
accept/reject correction to the drifted Gaussian proposal.

For each coordinate, it forms:

```julia
mu_f = R + D * dt * F(R)
mu_b = R' + D * dt * F(R')
```

and uses boundary-aware displacements to compute the forward and backward
Gaussian exponents. The resulting log acceptance ratio is:

```julia
log_ratio =
    -(sq_b - sq_f) / (4 * D * dt) +
    2 * (logpsi(R') - logpsi(R))
```

The move is accepted when:

```julia
log(rand()) < min(0.0, log_ratio)
```

If the move is rejected, the proposal buffers are overwritten with the original
state so that later code can treat accepted and rejected walkers uniformly.

If importance sampling is disabled, there is no Metropolis correction and the
proposal is treated as accepted immediately.

### 5. Update multiplicative weights

Each surviving walker multiplies its current weight by a short-time branching
factor:

```julia
logw = -0.5 * dt * ((Eold + Enew) - 2 * ET)
w_branch = exp(clamp(logw, -700.0, 700.0))
w_branch = min(w_branch, branch_cap)
```

and then:

```julia
weights[j] *= w_branch
```

Important details:

- the code uses the trapezoid average of old and new branch energies
- `ET` is the current reference energy
- the exponential is clamped for numerical safety
- the final branching factor is capped by `branch_cap`

For rejected guided moves, `Enew` has already been reset to `Eold`, so the
branch weight reduces to the expected diagonal form.

### 6. Swap buffers

At the end of the walker loop, the code swaps:

- `positions <-> proposal_positions`
- `state_data <-> proposal_data`

This makes the proposal buffers become the new current state without copying
all arrays back element by element.

### 7. Record step-level diagnostics

After the swap, `step!` pushes:

- mean walker weight: `sum(weights) / N`
- effective population:
  `(sum(weights)^2) / sum(weights.^2)`
- acceptance rate for that propagation step

Then it increments `sim.step`.

### 8. Reconfigure on schedule

If:

```julia
sim.step % sim.params.reconfiguration_interval == 0
```

the code resamples the weighted ensemble back to equal weights.

## Reconfiguration

GFMC keeps the number of walkers fixed by resampling indices from the current
weight distribution.

The helper `_normalized_weights(weights)` returns:

- `probs = weights / sum(weights)`
- `mean_weight = sum(weights) / N`

and throws an error if the total weight is not finite or is non-positive.

### Available reconfiguration schemes

- `MultinomialReconfiguration()`
- `SystematicReconfiguration()`

Both return a vector of source indices and the ensemble mean weight.

`_reconfigure!(sim)` then:

1. resamples source indices from `sim.weights`
2. copies the selected walker columns into `proposal_positions`
3. copies the selected batch-data entries into `proposal_data`
4. swaps the proposal buffers into the active state
5. resets every walker weight to `1.0`

After reconfiguration:

- the walker count is still exactly `targetN`
- the ensemble is now equally weighted again
- any information about the previous weight spread survives only through the
  already-recorded histories

`population_history` therefore stays constant even when the weighted ensemble
has effectively collapsed. To monitor collapse, use:

- `effective_population_history`
- `mean_weight_history`

## Reference-Energy Update

The reference energy update happens in `run_simulation!`, not inside `step!`.

It is applied only when both conditions are true:

- the current loop index `s` is greater than `nequil`
- the step just completed was a reconfiguration step

The update computes a block-averaged energy:

- if `energy_mean_history` is empty, use `estimate_energy(sim)`
- otherwise, average the last `energy_window` entries of
  `energy_mean_history`

Then it defines:

```julia
dt_block = dt * reconfiguration_interval
correction = (feedback / dt_block) * log(mean_weight)
ET = Eblock - correction
```

with the special case:

- if `feedback == 0`, the correction is skipped and `ET = Eblock`

Interpretation:

- if the recent average weight is larger than `1`, the feedback lowers `ET`
  so future branching factors shrink
- if the recent average weight is smaller than `1`, the feedback raises `ET`
  so future branching factors grow

This is the fixed-population analogue of population-control feedback in DMC.

## Run Loop

The public runner is:

```julia
run_gfmc!(sim; snapshot_steps=Int[], show_progress=false, progress_every=0,
          progress_label="", debug=false, debug_every=1, debug_io=stdout)
```

This is just an alias for `run_simulation!`.

### Initialization before the loop

Before any propagation step, the runner seeds the history arrays with the step
0 state:

- `mean_weight_history` gets `1.0`
- `effective_population_history` gets the current walker count
- `acceptance_history` gets `1.0`
- `record_state!` stores the initial energy, variance, population, and `ET`

### Main loop body

For each `s in 1:nsteps`, the runner:

1. calls `step!(sim)`
2. updates `ET` if equilibration is finished and this was a reconfiguration
   step
3. calls `record_state!(sim, s in snapshot_set)`
4. optionally emits progress output
5. optionally emits debug output

### What `snapshot_steps` actually controls

`snapshot_steps` only controls whether walker coordinates are stored in
`walker_positions_history`.

If `snapshot_steps` is left empty, positions are recorded at every step.

It does not control whether scalar histories are recorded. The following are
always recorded every step:

- `population_history`
- `energy_mean_history`
- `energy_variance_history`
- `ET_history`

If you pass `snapshot_steps=[0, params.nsteps]`, you still get full scalar
histories but only the initial and final coordinate snapshots.

## Energy Estimators and Recorded Histories

### `estimate_energy(sim)`

This computes the weighted mean of `sim.state_data.branch_energy`:

```julia
sum(weights[j] * E[j]) / sum(weights[j])
```

If the total weight is zero, it returns `0.0`.

### `estimate_energy_variance(sim)`

This computes the weighted variance around `estimate_energy(sim)`.

### `record_state!(sim, record_positions=true)`

Every call appends:

- current population
- current weighted mean branch energy
- current weighted branch-energy variance
- current `ET`

If `record_positions=true`, it also stores a deep-copied walker snapshot as a
`Vector{Vector{Float64}}`.

## Debug and Progress Reporting

The GFMC runner supports the same debug/progress pattern as the DMC and VMC
workflows.

Useful keywords:

- `show_progress=true`
- `progress_every=...`
- `progress_label="..."`
- `debug=true`
- `debug_every=...`
- `debug_io=...`

Debug lines include:

- `population`
- `ET`
- `mean_weight`
- `effective_population`
- `acceptance_rate`
- `energy_mean`
- `energy_var`
- whether the current step triggered reconfiguration

## Example Usage

### Generic Hamiltonian path

```julia
using Random
using System1D

bc = PeriodicBoundary1D(0.0, 1.0)
H = Hamiltonian(1, 0.5, R -> 0.5 * sum(abs2, R), bc)
trial = TrialWF(
    R -> -0.5 * sum(abs2, R),
    R -> -R,
    R -> -length(R),
)
guiding = ImportanceGuiding(trial, H)

params = GFMCParams(0.02, 100, 20, 256, 0.5, 0.2, 5, 5.0, 20)
walkers = [Walker([rand()]) for _ in 1:params.targetN]

sim = GFMCSim(H, params, walkers, MersenneTwister(7);
    guiding=guiding,
    nodepolicy=NoNode(),
    reconfiguration=SystematicReconfiguration(),
)

run_gfmc!(sim; snapshot_steps=[0, params.nsteps], debug=true)
```

### Generic Hamiltonian path with VMC initialization

```julia
using Random
using System1D

bc = PeriodicBoundary1D(0.0, 1.0)
H = Hamiltonian(1, 0.5, R -> 0.5 * sum(abs2, R), bc)
trial = TrialWF(
    R -> -0.5 * sum(abs2, R),
    R -> -R,
    R -> -length(R),
)

vmc_params = VMCParams(; dt=0.02, nsteps=100, targetN=256, ET0=0.0)
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

### Spinless fermion ring path

```julia
using Random
using System1D

model = SpinlessFermionRing1D(3, 1.0, 3.0, 1.0;
    lambda=-0.5,
    node_tol=1e-8,
    trig_eps=1e-10,
)
params = GFMCParams(0.003, 50, 10, 128, -0.2, 0.1, 2, 5.0, 10)
positions = sample_uniform_configurations(model, params.targetN, MersenneTwister(101))

sim = GFMCSim(model, params, positions, MersenneTwister(202);
    use_guiding=true,
    nodepolicy=FixedNode(),
)

run_gfmc!(sim; snapshot_steps=[0, params.nsteps])
```

### Spinless fermion ring path with VMC initialization

```julia
using Random
using System1D

model = SpinlessFermionRing1D(3, 1.0, 3.0, 1.0;
    lambda=-0.5,
    node_tol=1e-8,
    trig_eps=1e-10,
)
vmc_params = VMCParams(; dt=0.001, nsteps=50, targetN=128, ET0=0.0)
gfmc_params = GFMCParams(0.003, 50, 10, 128, -0.2, 0.1, 2, 5.0, 10)
seed_positions = sample_uniform_configurations(model, gfmc_params.targetN, MersenneTwister(101))

sim = run_gfmc_with_vmc_init(
    model,
    gfmc_params,
    seed_positions,
    vmc_params;
    vmc_rng=MersenneTwister(202),
    gfmc_rng=MersenneTwister(303),
    proposal=:gaussian,
    snapshot_steps=[0, gfmc_params.nsteps],
)
```

## Extending the GFMC Layer

To add a new physics model efficiently, the cleanest path is usually to define
a new `AbstractGFMCKernel` subtype and implement:

- `configuration_dimension(kernel)`
- `diffusion_constant(kernel)`
- `boundary(kernel)`
- `uses_importance_sampling(kernel)` if needed
- `uses_fixed_node(kernel)` if needed
- `evaluate_configuration_data!(kernel, data, X)`

The batch evaluator should fill `data` for every walker column in `X`.

If you do not need a specialized kernel, you can often stay with
`GenericGFMCKernel` by providing:

- a `Hamiltonian`
- a `TrialWF`
- `ImportanceGuiding(trial, H)` if guided sampling is desired
- a `FixedNode()` policy if sign-based killing is desired

To add a new resampling scheme, define a new subtype of
`AbstractGFMCReconfiguration` and add:

```julia
resample_indices(::MyScheme, rng, weights)
```

that returns `(indices, mean_weight)`.

The backend interface is also abstracted, but only `CPUBackend()` is currently
implemented.

## Current Limitations and Caveats

- The GFMC layer is narrower in scope than the DMC and VMC implementations.
- There is no generic high-level `run_gfmc(...)` convenience constructor/runner.
- The only built-in one-call convenience wrapper is `run_gfmc_with_vmc_init(...)`, which is specifically for VMC-backed warm starts.
- Only the CPU backend exists.
- Fixed-node enforcement only has real effect when the active kernel actually
  populates `sign`; the built-in unguided paths set `sign = 1.0`.
- Passing `Walker` objects into `GFMCSim(...)` discards their `phase` field,
  because GFMC stores walkers as a real-valued position matrix.
- GFMC uses coordinate wrapping only. It does not propagate the complex phase
  factors associated with twisted periodic boundaries, unlike the DMC walker
  path that uses `wrap_position_with_phase`.
- The recorded energy is the weighted mean of the branch-energy estimator, not
  a descendant-weight or forward-walking estimator.

## Reading the Code in Order

If you want to understand the implementation from top to bottom, read it in
this order:

1. `src/UseCases/gfmc/params.jl`
2. `src/UseCases/gfmc/abstract.jl`
3. `src/UseCases/gfmc/kernels.jl`
4. `src/UseCases/gfmc/sim.jl`
5. `src/UseCases/gfmc/reconfiguration.jl`
6. `src/UseCases/gfmc/run.jl`

That order mirrors the runtime dependency chain: parameterization, physics
adapter, state container, resampling policy, then the actual propagation loop.
