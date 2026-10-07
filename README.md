# Renormalization Group Flows of Spectral Dimensions on Hierarchical Lattices

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23201028.svg)](https://doi.org/10.5281/zenodo.23201028)

Code, data, figures and manuscript for

> **R. Chen**, *Renormalization Group Flows of Spectral Dimensions on Hierarchical Lattices: Critical Percolation,
> Resistance Exponents, and Transport Universality* (2026). DOI: [10.5281/zenodo.23201028](https://doi.org/10.5281/zenodo.23201028)

## Summary

- **Exact RG.** Bond percolation on the diamond hierarchical lattice (DHL) obeys p' = 2p² − p⁴. The fixed points
  0 (fragmented, d_s → 0) and 1 (intact lattice, d_s = 2) attract; p_c = (√5 − 1)/2 repels, with
  ν = ln 2 / ln(6 − 2√5) ≈ 1.6353. A non-integer spectral dimension appears only at a tuned critical point, never
  spontaneously.
- **Closed-form mass dimension.** For the critical spanning cluster, including dangling parts, the moments
  (E[M 1_c], E[M 1_not c]) obey an exact 2×2 recursion. Its Perron eigenvalue Λ = 7 − 2√5 + √(73 − 32√5) gives
  d_f = log₂ Λ = 1.89929.
- **Resistance exponent and spectral dimension.** The pole-to-pole resistance grows with Lyapunov exponent
  ζ = 0.8127 ± 0.0002, and the Einstein relation gives d_s = 2 d_f/(d_f + ζ) ≈ 1.401. Exact diagonalisation of
  generation-6 and generation-7 clusters, averaged over two log-periods, gives 1.41 and 1.44; the intact lattice
  gives 1.994 against the exact 2. So d_s ≈ 1.40–1.44, and the value 1.25 is excluded.
- **Crossovers and rewiring.** Off-critical clusters cross over to d_s = 0 or 2 beyond ξ ~ |p − p_c|^-ν. Newman–Watts
  shortcuts are relevant and raise d_s beyond τ ~ ξ² without producing a plateau.
- **Transport universality.** The DHL cluster has nearly the mass dimension of 2D percolation (1.8993 vs 91/48),
  but a different resistance exponent (0.813 vs 0.983) and spectral dimension (≈1.40 vs 1.317). The spectral
  dimension of "the" critical cluster also depends on the subgraph and Laplacian: the Hambly–Kumagai diffusion on
  the skeleton has 2θ/(θ+1) ≈ 1.68.

This is the seventh paper of a series:
[1](https://doi.org/10.5281/zenodo.23189575) holographic entanglement of fractal regions
([code](https://github.com/Ruqing1963/holographic-fractal-entanglement)),
[2](https://doi.org/10.5281/zenodo.23193512) multifractal large deviations in quantum entanglement
([code](https://github.com/Ruqing1963/multifractal-quantum-entanglement)),
[3](https://doi.org/10.5281/zenodo.23194419) uncomputable fractal dimensions and the Chaitin barrier
([code](https://github.com/Ruqing1963/chaitin-fractal-complexity)),
[4](https://doi.org/10.5281/zenodo.23195881) crossover from discrete to continuous quantum geometry
([code](https://github.com/Ruqing1963/quantum-geometry-crossover)),
[5](https://doi.org/10.5281/zenodo.23199122) thermodynamics of scale
([code](https://github.com/Ruqing1963/scale-thermodynamics)),
[6](https://doi.org/10.5281/zenodo.23199984) undecidable connectivity of computable fractal limit sets
([code](https://github.com/Ruqing1963/fractal-connectivity-undecidability)).

## Repository structure

| Path | Contents |
|---|---|
| `paper/` | LaTeX source and compiled PDF of the manuscript |
| `code/rg_spectral_dimension.py` | Single script producing every number and figure in the paper |
| `figures/` | Figure as vector PDF (used by the paper) and PNG |
| `data/` | Numerical data as CSV (metadata in `#` header lines) |
| `results/` | Console output of the script (the numbers quoted in the paper) |

## Reproducing the results

Requirements: Python ≥ 3.10 and the packages in `requirements.txt`
(tested with Python 3.12.4, NumPy 1.26.4, SciPy 1.13.1, Matplotlib 3.8.4). Runtime is about 15 minutes, dominated
by the dense diagonalisation of generation-7 clusters; random seeds are fixed.

```bash
pip install -r requirements.txt
python code/rg_spectral_dimension.py      # add --show to display the figure
```

Run from the repository root; the figure goes to `figures/` and data to `data/` (override with `RG_FIG_DIR`,
`RG_DATA_DIR`).

| Data file | Content | Paper |
|---|---|---|
| `rg_flow_trajectories.csv` | iterates of p' = 2p² − p⁴ from starting points near p_c | §2, Fig. 1 |
| `dhl_critical_cluster.csv` | crossing probability, mean mass (exact and sampled), resistance statistics per level | §3–4, Fig. 1 |
| `dhl_heat_kernel.csv` | d_s(τ) for intact, critical and off-critical DHL clusters | §4–5, Table 1, Fig. 1 |
| `newman_watts_ds.csv` | d_s(τ) for a ring with Newman–Watts shortcuts | §5, Fig. 1 |

To rebuild the paper (pdfLaTeX, two passes):

```bash
cd paper
pdflatex Chen_2026_RG_Spectral_Dimension.tex
pdflatex Chen_2026_RG_Spectral_Dimension.tex
```

## Citation

```bibtex
@misc{Chen2026RGSpectralDimension,
  author = {Chen, Ruqing},
  title  = {Renormalization Group Flows of Spectral Dimensions on Hierarchical Lattices: Critical Percolation,
            Resistance Exponents, and Transport Universality},
  year   = {2026},
  doi    = {10.5281/zenodo.23201028},
  url    = {https://doi.org/10.5281/zenodo.23201028}
}
```

## License

- **Code** (`code/`): [MIT License](LICENSE)
- **Manuscript, figures, data and results** (`paper/`, `figures/`, `data/`, `results/`):
  [CC BY 4.0](LICENSE-CC-BY-4.0.md)

## Contact

Ruqing Chen — GUT Geoservice Inc., Montreal — ruqing@hotmail.com
