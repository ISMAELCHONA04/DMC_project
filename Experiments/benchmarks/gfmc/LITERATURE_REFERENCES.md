# GFMC Literature References

This note records the primary-source papers used for the literature-backed GFMC
comparison layer in `generate_literature_report.jl`.

## Direct Analytical References

- Free particle ring and 1D harmonic oscillator:
  E. Schrodinger, *Quantisierung als Eigenwertproblem (Zweite Mitteilung)*,
  Annalen der Physik 79, 489-527 (1926).
  DOI: <https://doi.org/10.1002/andp.19263840602>
  Note:
  The ring benchmark uses the periodic rotor/free-wave sector as an inference
  from the same exact spectral framework. The harmonic-oscillator benchmarks
  use the exact `E0 = omega / 2` result directly.

- Two coupled harmonic oscillators:
  J. C. Vega, D. Ojeda-Guillen, and R. D. Mota,
  *Exact solution of the isotropic and anisotropic Hamiltonian of two coupled harmonic oscillators*,
  arXiv:2410.00021 (2024).
  Link: <https://arxiv.org/abs/2410.00021>
  Note:
  This supports the exact center-of-mass plus relative-mode decomposition used
  by the coupled-two-particle HO benchmark.

- One-dimensional hydrogen odd sector:
  R. Loudon, *One-Dimensional Hydrogen Atom*,
  American Journal of Physics 27, 649-655 (1959).
  DOI: <https://doi.org/10.1119/1.1934950>

- One-dimensional hydrogen revisited:
  G. Palma and U. Raff,
  *The one dimensional Hydrogen atom revisited*,
  Canadian Journal of Physics 84, 787-800 (2006).
  DOI: <https://doi.org/10.1139/P06-072>
  Note:
  The fixed-node hydrogen benchmark targets the odd-parity branch, so the
  comparison to `E = -1/2` is an explicit inference from the cited odd-state
  theory rather than a claim about the controversial singular even sector.

## Paper-Backed Exact Theory Evaluated At Benchmark Parameters

- Periodic potentials and Bloch theory:
  F. Bloch, *Uber die Quantenmechanik der Elektronen in Kristallgittern*,
  Zeitschrift fur Physik 52, 555-600 (1929).
  DOI: <https://doi.org/10.1007/BF01339455>

- Periodic-potential wave functions:
  J. C. Slater, *Wave Functions in a Periodic Potential*,
  Physical Review 51, 846-851 (1937).
  DOI: <https://doi.org/10.1103/PhysRev.51.846>

- Sinusoidal-potential Mathieu mapping:
  J. Yin and S. C. Erwin,
  *Noninteracting Electrons in a Prototypical One-Dimensional Sinusoidal Potential*,
  arXiv:2003.06647 (2020).
  Link: <https://arxiv.org/abs/2003.06647>
  Note:
  The cosine-lattice benchmark uses the exact Mathieu/Bloch theory established
  by these references. The numerical reference value in the report is the
  repository's evaluation of that exact theory for the benchmark parameter set.

## Theory Context For Fermionized And Hard-Core Limits

- Slater determinant construction:
  J. C. Slater, *The Theory of Complex Spectra*,
  Physical Review 34, 1293-1322 (1929).
  DOI: <https://doi.org/10.1103/PhysRev.34.1293>

- Tonks-Girardeau mapping:
  M. Girardeau, *Relationship between Systems of Impenetrable Bosons and Fermions in One Dimension*,
  Journal of Mathematical Physics 1, 516-523 (1960).
  DOI: <https://doi.org/10.1063/1.1703687>

- Classical hard-rod limit:
  L. Tonks, *The Complete Equation of State of One, Two and Three-Dimensional Gases of Hard Elastic Spheres*,
  Physical Review 50, 955-963 (1936).
  DOI: <https://doi.org/10.1103/PhysRev.50.955>

- Repulsive one-dimensional Bose gas:
  E. H. Lieb and W. Liniger,
  *Exact Analysis of an Interacting Bose Gas. I. The General Solution and the Ground State*,
  Physical Review 130, 1605-1616 (1963).
  DOI: <https://doi.org/10.1103/PhysRev.130.1605>

These papers are used as context for the bosonic TG-style and hard-core stress
benchmarks, but those cases are not plotted against a single published numeric
value because the repository benchmarks use custom external potentials and
smooth interactions rather than the exact paper parameter sets.
