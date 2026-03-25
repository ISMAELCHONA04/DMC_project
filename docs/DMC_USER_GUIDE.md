# DMC User Guide

## Overview

The DMC implementation lives under `src/UseCases/dmc/` and supports:

- plain diffusion Monte Carlo (`NoGuiding()`)
- importance-sampled DMC (`ImportanceGuiding`)
- optional fixed-node rejection (`FixedNode()`)
- periodic or open boundaries
- twisted periodic boundaries through a phase factor accumulated on each walker

`population_control_gain` controls the strength of the population-feedback term
in the DMC growth estimator. It defaults to `1.0`. The legacy keyword
`pop_gain` is still accepted for compatibility, but `population_control_gain`
is the public parameter name.

For the companion variational workflow, see `docs/VMC_USER_GUIDE.md`.

## Main Types

### `DMCParams`

Use either the positional constructor

```julia
DMCParams(dt, nsteps, nequil, targetN, ET0, population_control_gain, branch_cap, nblocks)
```

or the keyword constructor

```julia
DMCParams(; dt, nsteps, targetN, nequil=0, ET0=0.0, population_control_gain=1.0, branch_cap=8, nblocks=50)
```

Parameters:

- `dt`: imaginary-time step size
- `nsteps`: number of DMC propagation steps
- `nequil`: number of equilibration steps before `update_ET!` starts adjusting the reference energy
- `targetN`: target walker population used by the growth estimator
- `ET0`: initial reference-energy guess
- `population_control_gain`: feedback strength multiplying the population-control correction in `update_ET!`
- `branch_cap`: hard cap on the branching factor for each walker in one step
- `nblocks`: averaging window used when smoothing `ET_history` inside `update_ET!`

Validation rules:

- `dt > 0`
- `nsteps >= 1`
- `0 <= nequil <= nsteps`
- `targetN >= 1`
- `population_control_gain >= 0`
- `branch_cap >= 1`
- `nblocks >= 1`

### `Walker`

Each walker stores:

- `position(w)`: the wrapped configuration vector
- `phase(w)`: a complex phase factor, defaulting to `1 + 0im`

The phase is updated automatically when a move crosses a twisted periodic
boundary.

### `PeriodicBoundary1D`

Periodic boundaries now accept an optional twist:

```julia
bc = PeriodicBoundary1D(0.0, L; twist=pi / 3)
```

or a per-coordinate twist vector:

```julia
bc = PeriodicBoundary1D(0.0, L; twist=[pi / 3, pi / 3, pi / 3])
```

Wrapping utilities:

- `wrap_position(bc, R)` returns only the wrapped position
- `wrap_position_with_phase(bc, R)` returns `(wrapped_R, phase_factor)`

For open boundaries, or periodic boundaries with `twist=0.0`, the phase factor
is always `1 + 0im`.

## Defining the Physics Input

This section covers the four objects that define a DMC problem:

- the Hamiltonian
- the boundary condition
- the trial wave function, if you use importance sampling
- the nodal policy, if you want fixed-node DMC

### Defining a Hamiltonian

The base Hamiltonian constructor is

```julia
Hamiltonian(N, D, V, bc=OpenBoundary1D())
```

with:

- `N`: number of coordinates in the configuration vector `R`
- `D`: diffusion constant
- `V`: potential callback with signature `R -> scalar`
- `bc`: boundary-condition object

The potential function always receives the full walker configuration `R` as a
vector. For a one-particle problem, `R` is still a one-element vector. For a
many-particle problem, `R[i]` is the `i`th coordinate.

Simple example:

```julia
H = Hamiltonian(1, 0.5, R -> 0.5 * R[1]^2)
```

Two-particle example:

```julia
H = Hamiltonian(2, 0.5, R -> 0.5 * (R[1]^2 + R[2]^2) + 0.1 * (R[1] - R[2])^2)
```

For periodic systems, define the potential in a way that is consistent with the
cell geometry. If your potential depends on separations, use the boundary
helpers such as `displacement(bc, x1, x2)` or `minimum_image(bc, dx)` rather
than raw differences.

Example with a ring:

```julia
bc = PeriodicBoundary1D(0.0, 1.0)
H = Hamiltonian(2, 0.5, R -> begin
    dx = displacement(bc, R[1], R[2])
    0.25 * dx^2
end, bc)
```

### Choosing boundary conditions

The two boundary types are:

- `OpenBoundary1D()`
- `PeriodicBoundary1D(xmin, xmax; twist=0.0)`

Use open boundaries when the coordinate is not wrapped. Use periodic
boundaries when the simulation cell is a ring or finite periodic interval.

Examples:

```julia
bc_open = OpenBoundary1D()
bc_periodic = PeriodicBoundary1D(0.0, 10.0)
bc_twisted = PeriodicBoundary1D(0.0, 10.0; twist=pi / 3)
```

Notes:

- periodic wrapping is applied coordinate-by-coordinate
- if `twist` is a scalar, every coordinate uses the same twist angle
- if `twist` is a vector, its length must match the configuration dimension
- the twist phase is accumulated when a walker crosses the boundary, but the
  potential callback itself still just sees a coordinate vector

Useful helpers when defining periodic models:

- `wrap_coordinate(bc, x)`
- `wrap_position(bc, R)`
- `wrap_position_with_phase(bc, R)`
- `minimum_image(bc, dx)`
- `displacement(bc, x1, x2)` or `displacement(bc, R1, R2)`

### Defining a trial wave function

If you use importance sampling, define a `TrialWF` and wrap it in
`ImportanceGuiding(trial, H)`.

The constructor is:

```julia
TrialWF(logpsi, gradlogpsi, lapllogpsi)
TrialWF(logpsi, gradlogpsi, lapllogpsi, signpsi)
```

Expected callbacks:

- `logpsi(R)`: returns `log(abs(psi_T(R)))`
- `gradlogpsi(R)`: returns the gradient of `log(abs(psi_T))` as a vector with the same length as `R`
- `lapllogpsi(R)`: returns the Laplacian of `log(abs(psi_T))`
- `signpsi(R)`: returns `-1.0`, `0.0`, or `1.0`

If you omit `signpsi`, it defaults to `+1` everywhere. That is appropriate for
bosonic or nodeless trial states, but not for fixed-node calculations.

The DMC code uses the trial wave function as follows:

- drift term: `2 * gradlogpsi(R)`
- local energy: `-D * (lapllogpsi(R) + sum(abs2, gradlogpsi(R))) + V(R)`
- node handling: `signpsi(R)`

Gaussian example for a 1D harmonic oscillator:

```julia
trial = TrialWF(
    R -> -0.5 * R[1]^2,
    R -> [-R[1]],
    R -> -1.0,
)

guiding = ImportanceGuiding(trial, H)
```

Two-coordinate example:

```julia
trial = TrialWF(
    R -> -0.25 * (R[1]^2 + R[2]^2),
    R -> [-0.5 * R[1], -0.5 * R[2]],
    R -> -1.0,
)
```

### Defining nodal structure

Nodal structure only matters when you use `FixedNode()`. The available policies
are:

- `NoNode()`: ignore the sign structure entirely
- `FixedNode()`: reject moves that cross a node

The fixed-node rule in this code is:

- evaluate `signpsi` at the old and proposed positions
- reject the move if the signs differ
- also reject the move if either sign is exactly `0.0`

That means your `signpsi` callback should return `0.0` in a small region around
the node if you want near-node proposals to be rejected cleanly.

Odd-parity 1D example with a node at `x = 0`:

```julia
trial = TrialWF(
    R -> log(abs(R[1])) - 0.5 * R[1]^2,
    R -> [inv(R[1]) - R[1]],
    R -> -(inv(R[1])^2) - 1.0,
    R -> sign(R[1]),
)

guiding = ImportanceGuiding(trial, H)
nodepolicy = FixedNode()
```

In practice you usually do not want the exact singular form above at `x = 0`.
Instead, define `signpsi` with a small tolerance and keep `logpsi`,
`gradlogpsi`, and `lapllogpsi` numerically safe near the node.

Example with a tolerance:

```julia
node_tol = 1e-8

trial = TrialWF(
    R -> log(max(abs(R[1]), node_tol)) - 0.5 * R[1]^2,
    R -> abs(R[1]) < node_tol ? [0.0] : [inv(R[1]) - R[1]],
    R -> abs(R[1]) < node_tol ? -1.0 : (-(inv(R[1])^2) - 1.0),
    R -> abs(R[1]) < node_tol ? 0.0 : sign(R[1]),
)
```

Practical rules:

- use `NoGuiding()` if you do not want a trial state at all
- use `ImportanceGuiding(trial, H)` if you want drift and local-energy
  branching from a trial state
- use `FixedNode()` only when the sign structure is physically meaningful and
  your `signpsi` callback is well defined
- if you use `FixedNode()` together with `NoGuiding()`, there is effectively no
  nodal structure because `signpsi(::NoGuiding, R)` is always `1.0`

## Running a Simulation

### Recommended one-call interface

```julia
using Random
using System1D

bc = PeriodicBoundary1D(-5.0, 5.0; twist=pi / 4)
H = Hamiltonian(1, 0.5, R -> 0.5 * R[1]^2, bc)
params = DMCParams(; dt=0.01, nsteps=200, nequil=20, targetN=256, ET0=0.5, population_control_gain=1.0, branch_cap=4, nblocks=20)
positions = [[randn()] for _ in 1:params.targetN]

sim = run_dmc(H, params, positions;
    rng=MersenneTwister(42),
    snapshot_steps=[0, params.nsteps],
    show_progress=true,
)
```

`run_dmc(...)` accepts:

- a vector of raw position vectors, such as `[[x1], [x2], ...]`
- a single raw position vector, which becomes a one-walker simulation
- a matrix whose columns are walker configurations
- an explicit vector of `Walker` objects

Keyword arguments:

- `rng`
- `guiding`
- `nodepolicy`
- `step`
- `snapshot_steps`
- `show_progress`
- `progress_every`
- `progress_label`
- `debug`
- `debug_every`
- `debug_io`

### Low-level interface

If you want to build the simulation object first and run later:

```julia
sim = DMCSim(H, params, positions, MersenneTwister(42); guiding=NoGuiding(), nodepolicy=NoNode())
run_simulation!(sim; snapshot_steps=[0, 50, 100])
```

This is still the core interface used internally. The new `run_dmc(...)` helper
just wraps this flow.

## Important Functions

### `DMCSim(...)`

Constructs the simulation state. Initial positions are wrapped immediately. If
the boundary is twisted and an initial configuration lies outside the cell, the
walker's stored phase is updated during that initial wrap.

### `initialize!(sim, configs)`

Resets an existing `DMCSim` with new walkers or raw positions and clears all
histories.

### `step!(sim)`

Runs one propagation step:

- proposes a move
- wraps it back into the cell if needed
- accumulates any twist phase on the walker
- applies node rejection if enabled
- computes the branching factor
- replicates or deletes walkers

### `run_simulation!(sim; ...)`

Runs `sim.params.nsteps` iterations, records diagnostics, and updates `ET`
after equilibration.

Diagnostics recorded in the simulation object:

- `ET_history`
- `population_history`
- `energy_mean_history`
- `energy_variance_history`
- `walker_positions_history`

Debug mode:

- set `debug=true` to print a per-step line with timing, population, `ET`, and
  current energy statistics
- use `debug_every` to thin the output if you only want every `n`th step
- use `debug_io` if you want to capture the diagnostics in an `IOBuffer` or a file

### `update_ET!(sim)`

Applies the population-control update

```julia
ET <- Eblock - (population_control_gain / dt) * log(Nt / N0)
```

where `Eblock` is the mean of the last `nblocks` entries in `ET_history`.

### `propose_move(sim, walker, guiding)`

Returns only the wrapped trial position. This remains useful for inspection and
testing. The internal DMC step uses an extended move object that also carries
the twist phase factor.

## Twisted Boundary Conditions

With a twisted periodic boundary, every time a coordinate crosses the periodic
cell by one full winding, the walker picks up a factor

```julia
exp(im * theta)
```

for positive winding, or the inverse factor for negative winding.

For a configuration with several coordinates, the total phase is the product of
the phase contributions from each wrapped coordinate.

Example:

```julia
bc = PeriodicBoundary1D(0.0, 1.0; twist=pi / 2)
wrapped, phase_factor = wrap_position_with_phase(bc, [1.25])

# wrapped == [0.25]
# phase_factor == cis(pi / 2)
```

Important scope note:

- the current implementation tracks the twist phase on each walker and carries
  it through branching
- the existing DMC drift, local-energy, and branching formulas remain the same
  real-valued formulas as before
- if you need a fully phase-coupled or fixed-phase DMC formulation, that would
  require extending the Hamiltonian and guiding-wavefunction interfaces beyond
  the current real-valued API

## Typical Workflow

1. Define a Hamiltonian and boundary condition.
2. Build `DMCParams`.
3. Prepare initial walker positions.
4. Run with `run_dmc(...)` or `DMCSim(...)` plus `run_simulation!(...)`.
5. Inspect `energy_mean_history`, `ET_history`, `population_history`, and the
   final walker phases through `phase.(sim.walkers)`.
