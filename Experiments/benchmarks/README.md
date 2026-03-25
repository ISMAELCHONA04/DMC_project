# Benchmarks

This directory contains reproducible benchmark harnesses built on top of the
public experiment systems.

Benchmarks are intended to answer two separate questions:

- Does a code path still run end-to-end after changes?
- Does it still reproduce the expected physics and numerical behavior?

The current curated benchmark surface is:

- [`gfmc/`](./gfmc/): fixed-population GFMC smoke tests, parameter sweeps, and
  recommended accuracy runs across the published GFMC systems.
