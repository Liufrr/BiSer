# BiSer

Benchmarking code and analysis pipeline for **BiSer** (Bidirectional Seriation), a discretization-aware joint seriation method for genomic matrices. This repository accompanies the manuscript and provides all code needed to reproduce the simulation studies and real-data analyses.

## Overview

BiSer is an algorithm that jointly orders rows and columns of a non-negative matrix by:

1. **SVD-based bipartite embedding**: Rows and columns are jointly embedded into a shared spectral space via normalized SVD of the bipartite adjacency matrix.
2. **TSP-based path optimization**: A Traveling Salesman Problem (TSP) solver is applied to the similarity matrix of the embedded coordinates to find an optimal linear ordering.
3. **Boundary detection**: Structural boundaries are automatically identified from the diagonal profile of the reordered similarity matrix using peak/valley detection (via `scipy.signal.find_peaks`).

## Repository Structure

```
BiSer/
├── R/                                  # Core R functions
│   ├── methods.R                       # BiSer and all comparison methods
│   ├── data_generation.R               # Simulation data generators
│   ├── metrics.R                       # Evaluation metrics (ARI, NMI, etc.)
│   ├── visualization.R                 # Shared plotting utilities
│   ├── feature_selection.R             # ANOVA-based gene selection (supervised)
│   └── feature_selection_unsupervised.R  # HVG-based gene selection (unsupervised)
│
├── python/                             # Python utilities
│   ├── auto_boundaries.py              # Boundary detection via scipy
│   └── requirements.txt                # Python dependencies
│
├── simulations/                        # Simulation benchmarking
│   ├── run_simulation.R                # Unified simulation runner
│   ├── config_low_noise.R              # Low noise parameter configs
│   ├── config_mid_noise.R              # Mid noise parameter configs
│   ├── config_high_noise.R             # High noise parameter configs
│   └── sensitivity_analysis.R          # Hyperparameter sensitivity analysis
│
├── real_data/                          # Real dataset analyses
│   ├── khan_analysis.R                 # Khan SRBCT dataset analysis
│   ├── li_analysis.R                   # Li single-cell dataset analysis
│   ├── nci60_analysis.R                # NCI60 cancer cell line analysis
│   ├── analysis_yan_unsupervised.R     # Yan embryo dataset analysis
│   ├── compare_metrics.R              # Cross-dataset metric comparison
│   └── select_markers.R               # Candidate gene marker selection
│
├── scalability/                        # Computational complexity
│   └── benchmark_scalability.R         # Timing experiments
│
├── figures/                            # Figure generation scripts
│   ├── fig2_simulation_heatmap.R       # Figure 2: Simulation heatmaps
│   ├── fig3_simulation_ridge.R         # Figure 3: Ridge distributions
│   ├── fig4_discretization_penalty.R   # Figure 4: Discretization penalty
│   ├── fig5_boxplot_main.R             # Figure 5: 4-scenario boxplots
│   ├── fig6_yan_combined.R             # Figure 6: Yan combined (Panel B)
│   └── fig7_scalability.R              # Figure 7: Scalability benchmarks
│
├── data/                               # Data directory (see data/README.md)
├── output/                             # Generated figures and results
├── README.md
├── LICENSE
├── .gitignore
└── .gitattributes
```

## Requirements

### R (>= 4.0)

Install all required R packages:

```r
# CRAN packages
install.packages(c(
  # Seriation and clustering
  "TSP", "seriation", "mclust", "biclust", "NMF", "clue",

  # Visualization
  "ggplot2", "ggridges", "ggrepel", "circlize",
  "viridis", "cowplot", "patchwork", "gridExtra",

  # Data manipulation
  "dplyr", "tidyr", "tibble",

  # Metrics
  "aricode",

  # Python interface
  "reticulate",

  # Other
  "Matrix", "abind", "gtools", "igraph"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "scran"))
```

### Python (>= 3.8)

```bash
pip install -r python/requirements.txt
# or: pip install numpy scipy
```

### Python configuration for reticulate

Make sure the R `reticulate` package can find your Python installation:

```r
library(reticulate)
# Option 1: Use a specific Python
use_python("/path/to/python")

# Option 2: Use a conda environment
use_condaenv("your_env_name")
```

## Usage

**Important**: All scripts assume the working directory is set to the repository root.

```r
setwd("path/to/BiSer")
```

### Running simulations

```r
# Load core functions
source("R/methods.R")
source("R/data_generation.R")
source("R/metrics.R")

# Run low-noise simulations
source("simulations/config_low_noise.R")
source("simulations/run_simulation.R")
```

### Generating figures

```r
# After simulation data is available in output/
source("figures/fig2_simulation_heatmap.R")
source("figures/fig3_simulation_ridge.R")
source("figures/fig5_boxplot_main.R")
```

### Running sensitivity analysis

```r
# Requires saved simulation results (.RData files from run_simulation.R)
source("simulations/sensitivity_analysis.R")
# Outputs: output/FigS_sensitivity.png, output/sensitivity_raw_results.csv
```

### Running real-data analysis

```r
# Ensure data files are in data/ (see data/README.md)
source("real_data/khan_analysis.R")
source("real_data/li_analysis.R")
source("real_data/nci60_analysis.R")
source("real_data/analysis_yan_unsupervised.R")
```

## Benchmarked Methods

| Method | Type | Implementation |
|--------|------|----------------|
| **BiSer** | Joint seriation | Custom (this study) |
| BiSer_SVD | Ablation (no path optimization) | Custom (this study) |
| BiSer_TSP | Ablation (no joint embedding) | Custom (this study) |
| Spectral | Independent seriation | `seriation::seriate(method="spectral")` |
| Heatmap (OLO) | Independent seriation | `seriation::seriate(method="OLO")` |
| MESBC | Biclustering | Custom (Liu et al. 2024) |
| NMF | Biclustering | `NMF::nmf` (Brunet et al. 2004) |

## Evaluation Metrics

- **ARI** (Adjusted Rand Index)
- **NMI** (Normalized Mutual Information)
- **Purity**
- **Discretization penalty** (ARI loss from continuous order to discrete labels)
- **Spearman correlation** (ordering quality for developmental data)

## Data

See [`data/README.md`](data/README.md) for instructions on obtaining the required datasets.


## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
