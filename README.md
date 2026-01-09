# DMC_project

Project layout:

- `Project.toml`, `Manifest.toml`
- `src/` - package source
  - `DMC.jl` - entry module
  - `Domain/` - physics types
  - `UseCases/` - algorithms (propagate, branch, estimators)
  - `Infrastructure/` - RNG, IO, parallel, logging
  - `Interfaces/` - CLI / config / adapters
- `test/` - tests (contains `runtests.jl`)
- `scripts/` - runnable experiments (HO, H atom, etc.)

Quick test (from project root):

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Usage examples (plain vs importance sampling):

```julia
using Random
using System1D

H = Hamiltonian(1, 0.5, R -> 0.5 * R[1]^2)
params = DMCParams(0.01, 100, 10, 200, 0.5, 0.1, 2, 10)
walkers_plain = [Walker([randn()]) for _ in 1:params.targetN]
walkers_guided = [Walker([randn()]) for _ in 1:params.targetN]
rng = MersenneTwister(42)

sim_plain = DMCSim(H, params, walkers_plain, System1D.NoGuiding(), rng)

trial = System1D.TrialWF(
    R -> -0.5 * R[1]^2,
    R -> [-R[1]],
    R -> -1.0
)
guiding = System1D.ImportanceGuiding(trial, H)
sim_guided = DMCSim(H, params, walkers_guided, guiding, rng)
```
