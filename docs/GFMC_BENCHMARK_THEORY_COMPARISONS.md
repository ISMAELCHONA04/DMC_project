# GFMC Benchmark Theory Comparisons

This document collects the benchmark cases whose outputs can be compared directly to paper-backed theory references already cited in the repository.

Regenerate the benchmark outputs and these curated figures from the repository root with:

```bash
julia --project=. Experiments/benchmarks/gfmc/generate_theory_comparison_gallery.jl
```

## Curated Aggregate Figures

- [Final energies vs theory](assets/gfmc_theory_comparisons/final_energy_vs_theory.png)
- [Final energy error vs theory](assets/gfmc_theory_comparisons/final_energy_error_vs_theory.png)

## Included Cases

| Case | History/Error Plot | Density Plot | Theory Reference |
|---|---|---|---|
| Free Ring | [history/error](assets/gfmc_theory_comparisons/free_ring_history.png) | [density](assets/gfmc_theory_comparisons/free_ring_density.png) | [Schrodinger 1926 (II)](https://doi.org/10.1002/andp.19263840602) |
| HO Unguided | [history/error](assets/gfmc_theory_comparisons/ho_unguided_history.png) | n/a | [Schrodinger 1926 (II)](https://doi.org/10.1002/andp.19263840602) |
| HO Guided | [history/error](assets/gfmc_theory_comparisons/ho_guided_history.png) | n/a | [Schrodinger 1926 (II)](https://doi.org/10.1002/andp.19263840602) |
| Coupled 2-HO Guided | [history/error](assets/gfmc_theory_comparisons/coupled_2ho_guided_history.png) | n/a | [Vega et al. 2024](https://arxiv.org/abs/2410.00021) |
| Odd HO Fixed Node | [history/error](assets/gfmc_theory_comparisons/odd_ho_fixed_node_history.png) | n/a | [Schrodinger 1926 (II)](https://doi.org/10.1002/andp.19263840602) |
| 1D Hydrogen Fixed Node | [history/error](assets/gfmc_theory_comparisons/hydrogen_fixed_node_history.png) | n/a | [Loudon 1959](https://doi.org/10.1119/1.1934950), [Palma and Raff 2006](https://doi.org/10.1139/P06-072) |
| 1D Hydrogen Unguided | [history/error](assets/gfmc_theory_comparisons/hydrogen_unguided_history.png) | n/a | [Loudon 1959](https://doi.org/10.1119/1.1934950), [Palma and Raff 2006](https://doi.org/10.1139/P06-072) |
| Cosine Lattice Unguided | [history/error](assets/gfmc_theory_comparisons/cosine_lattice_unguided_history.png) | [density](assets/gfmc_theory_comparisons/cosine_lattice_unguided_density.png) | [Bloch 1929](https://doi.org/10.1007/BF01339455), [Slater 1937](https://doi.org/10.1103/PhysRev.51.846), [Yin and Erwin 2020](https://arxiv.org/abs/2003.06647) |
| Cosine Lattice Guided | [history/error](assets/gfmc_theory_comparisons/cosine_lattice_guided_dt0p0010_history.png) | [density](assets/gfmc_theory_comparisons/cosine_lattice_guided_dt0p0010_density.png) | [Bloch 1929](https://doi.org/10.1007/BF01339455), [Slater 1937](https://doi.org/10.1103/PhysRev.51.846), [Yin and Erwin 2020](https://arxiv.org/abs/2003.06647) |

## Final-Value Summary

| Benchmark | Computed | Theory | Error |
|---|---:|---:|---:|
| HO Unguided | 0.46103526 +/- 1.61e-03 | 0.50000000 | -3.896e-02 |
| Cosine Lattice Unguided | -1.92148817 +/- 2.96e-03 | -1.90802366 | -1.346e-02 |
| Coupled 2-HO Guided | 1.27459667 +/- 2.87e-17 | 1.27459667 | -7.327e-15 |
| Odd HO Fixed Node | 1.50000000 +/- 9.39e-16 | 1.50000000 | -6.661e-16 |
| Free Ring | 0.00000000 +/- 0.00e+00 | 0.00000000 | 0.000e+00 |
| HO Guided | 0.50000000 +/- 0.00e+00 | 0.50000000 | 0.000e+00 |
| 1D Hydrogen Fixed Node | -0.50000000 +/- 3.81e-18 | -0.50000000 | 0.000e+00 |

## Notes

- Only cases with direct paper-backed theory references are included here. The periodic-ion soft-Coulomb benchmarks remain exact within the repository model, but not against a matching published parameter set.
- The cosine-lattice guided comparison uses the smallest time-step sweep run (`dt = 0.0010`) because the locked final benchmark for that case keeps only the unguided branch.
- The unguided hydrogen case is included as a stress comparison against the odd-state 1D hydrogen reference; it is expected to fail badly, and that mismatch is the point of the plot.
