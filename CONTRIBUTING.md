# Contributing

This repository is maintained as a small, readable research codebase rather
than a large framework. Contributions should keep that character intact.

## Setup

Instantiate the root environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the package tests before submitting changes:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

If you work on notebooks, you may also want the optional notebook environment:

```bash
julia --project=Experiments/env -e 'using Pkg; Pkg.instantiate()'
```

## Expectations

- Keep public APIs stable unless a change is necessary and documented.
- Update tests and the relevant user guide in the same patch as behavioral changes.
- Prefer small, explicit helpers over dense multi-purpose abstractions.
- Add comments and docstrings only where they explain intent, math, or a subtle implementation choice.

## Experiment policy

- Scripts are preferred for reproducible public entrypoints.
- Notebooks should be committed without saved outputs.
- Generated figures, animations, CSV tables, and similar byproducts belong under gitignored `outputs/` directories.
- Do not commit generated LaTeX build products, notebook figure dumps, or machine-specific files.

## Style notes

- The repository uses plain Julia source files with explicit includes instead of metaprogrammed registration.
- Boundary handling, guiding, and node policies are kept separate on purpose so method layers can compose them explicitly.
- When adding a new experiment, document it in `Experiments/systems/README.md`.
