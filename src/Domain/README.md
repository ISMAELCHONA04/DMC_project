# Domain Layer

The `src/Domain/` directory contains the method-agnostic physics layer used by
the DMC, VMC, and GFMC implementations.

## Main responsibilities

- boundary and geometry policies (`OpenBoundary1D`, `PeriodicBoundary1D`)
- Hamiltonian containers and accessors
- walker state representation
- trial-wavefunction callbacks for guiding and fixed-node logic
- boundary-crossing phase handling for twisted periodic systems
- built-in model definitions such as `SpinlessFermionRing1D`

## Design intent

Code in this layer should stay as independent as possible from any particular
Monte Carlo algorithm. The DMC, VMC, and GFMC implementations consume these
types, but the domain layer itself should not depend on method-specific state
machines or history recording.

## Files

- `boundaries1D.jl`: geometry helpers, minimum-image logic, and twist phases
- `hamiltonian.jl`: `Hamiltonian` container and lightweight accessors
- `walker.jl`: `Walker` state type plus convenience constructors/accessors
- `trialwf.jl`: `TrialWF` callback container
- `guiding.jl`: `NoGuiding`, `ImportanceGuiding`, `drift`, and `local_energy`
- `nodepolicy.jl`: `NoNode` and `FixedNode`
- `spinless_fermion_ring_1d.jl`: specialized periodic-ring model used by the fermion workflows
