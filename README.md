# BiSer

This repository contains the implementation and analysis code corresponding to the BiSer0608 submitted manuscript. BiSer jointly orders the row and column nodes of a non-negative genomic matrix by normalized bipartite SVD, a shared spectral embedding, the joint Gram similarity matrix `K = Z %*% t(Z)`, and global path optimization. Discrete labels are induced from a Gaussian-smoothed diagonal profile of the reordered joint similarity matrix.

## Submitted simulation design

The primary benchmark is defined in `simulations/config_submission.R`. It contains the complete factorial design reported in Supplementary Tables S1–S3:

- Gaussian, Poisson, and negative binomial data-generating models
- mutually exclusive and overlapping biclusters
- low and high noise
- 250×150 and 1000×650 matrices
- 100 independent replicates in each of 24 conditions

The small overlapping design uses row overlaps `(0, 0, 8, 3, 20)` and column overlaps `(0, 0, 10, 20, 15)`. The large design uses row overlaps `(0, 10, 20, 30, 20)` and column overlaps `(0, 20, 10, 20, 15)`.

Run the complete submitted benchmark from the repository root:

```r
source("simulations/run_submission.R")
```

The low- and high-noise subsets can be run separately:

```r
source("R/methods.R")
source("R/metrics.R")
source("simulations/config_low_noise.R")  # or config_high_noise.R
source("simulations/run_simulation.R")
```

`config_mid_noise.R` is retained only as a legacy exploratory design and is not part of the submitted 24-condition benchmark.

## Methods

The submitted comparison contains exactly seven methods:

| Stored name | Manuscript name | Definition |
|---|---|---|
| `biser` | BiSer | Shared nontrivial SVD embedding and joint TSP path |
| `bs` | BiSer_SVD | Independent ordering by the leading nontrivial singular vectors |
| `tsp_seri` | BiSer_TSP | Independent row and column TSP paths in correlation-distance spaces |
| `spec_seri` | Spectral | Independent spectral seriation |
| `Heatmap` | Heatmap | Independent optimal leaf ordering (OLO) |
| `MESBC` | MESBC | Spectral biclustering with modularity selection over K=2,…,10 |
| `NMF` | NMF | Brunet NMF with modularity selection over K=2,…,10 |

BiSer removes zero-degree nodes before normalization, excludes the degree-related trivial singular component, retains all remaining singular components, and solves the joint path using 100 repeated nearest-neighbor starts followed by two-opt refinement.

## Boundary detection and evaluation

The boundary profile averages a symmetric diagonal band while excluding self-similarities. The profile is smoothed with a one-dimensional Gaussian filter (`sigma=3` by default), after which valleys are retained using minimum prominence and spacing constraints.

Separate row and column scores are combined by an arithmetic mean weighted by the corresponding numbers of row and column nodes. The saved metric matrices contain ARI, NMI, purity, accuracy, precision, recall, F1, recovery, and relevance, together with:

- `ARI_pre`: oracle segmentation of the recovered continuous order using the true number and sizes of groups
- `ARI_post`: ARI after data-driven boundary detection
- `discretization_penalty`: `ARI_pre - ARI_post`

The discretization penalty is defined for BiSer and the four ordering or ablation methods. MESBC and NMF return discrete biclusters directly, so their boundary-induced penalty is left undefined.

The sensitivity analysis evaluates the submitted 175-combination grid while holding Gaussian smoothing fixed:

```r
source("simulations/sensitivity_analysis.R")
```

## Figures

After generating the required simulation outputs:

```r
source("figures/fig2_simulation_heatmap.R")
source("figures/fig3_simulation_ridge.R")
source("figures/fig4_discretization_penalty.R")
source("figures/fig5_boxplot_main.R")
```

Figures 2 and 3 select the replicate nearest to the median BiSer ARI and use the exact disturbed matrix saved by the runner. Figure 4 plots pre-boundary versus post-boundary ARI, rather than a low-noise versus high-noise contrast.

## Requirements

R packages: `TSP`, `seriation`, `scran`, `igraph`, `biclust`, `NMF`, `mclust`, `aricode`, `clue`, `reticulate`, and the plotting packages used by individual figure scripts.

Python packages:

```bash
pip install -r python/requirements.txt
```

All scripts assume that the working directory is the repository root.
