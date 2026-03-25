# Documentation

This directory contains the user-facing and developer-facing documentation for
the repository.

## Guides

- `DMC_USER_GUIDE.md`: public DMC API, physics inputs, and run-loop behavior
- `VMC_USER_GUIDE.md`: VMC API, proposal modes, and estimator conventions
- `GFMC_USER_GUIDE.md`: fixed-population GFMC API, kernels, and warm starts
- `GFMC_BENCHMARK_THEORY_COMPARISONS.md`: curated simulation-vs-theory figures for paper-backed GFMC benchmarks
- `GFMC_CODE_WALKTHROUGH.tex`: LaTeX source for a deeper implementation walkthrough

## Building the LaTeX walkthrough

The repository keeps the `.tex` source, but not generated LaTeX byproducts.
From the repository root:

```bash
mkdir -p docs/build
latexmk -pdf -output-directory=docs/build docs/GFMC_CODE_WALKTHROUGH.tex
```

Generated files will be written under `docs/build/` and should not be
committed.

## Related references

- top-level project overview: [`../README.md`](../README.md)
- experiment inventory: [`../Experiments/README.md`](../Experiments/README.md)
- collaborator workflow: [`../CONTRIBUTING.md`](../CONTRIBUTING.md)
