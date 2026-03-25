# Use-Case Layer

The `src/UseCases/` directory contains the method-specific runtime logic built
on top of the shared domain layer.

## Subdirectories

- `common/`: helpers reused by multiple methods, including debug logging,
  energy summaries, snapshot recording, and VMC-backed warm starts
- `dmc/`: diffusion Monte Carlo state, stepping logic, branching, and plotting
- `vmc/`: variational Monte Carlo state, proposal mechanics, and estimators
- `gfmc/`: fixed-population GFMC kernels, reconfiguration, state, and run loop

## Design intent

Changes in this layer should preserve the public API described in the method
guides under `docs/`. If a change affects user-visible behavior, the
corresponding guide and tests should be updated in the same patch.
