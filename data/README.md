# Data Directory

This directory should contain the data files required to run the real data analyses.
Data files are not included in the repository due to size constraints.

## Expected Data Files

### Simulation Data

Simulation data is generated automatically by the scripts in `simulations/`.
Generated `.RData` files will be saved to `output/`.

### Real Datasets

The following `.RData` files are expected for real data analysis.
Each file should contain the named objects listed below.

#### `Khan.RData`
- **Source**: Khan et al. (2001), Small Round Blue Cell Tumors (SRBCT)
- **Reference**: Khan J, et al. *Nature Medicine* 2001;7:673-679
- **Objects**: Khan dataset from the `ISLR` R package
- **Obtaining**: Available via `ISLR::Khan` or from the original publication

#### `li.RData`
- **Source**: Li et al. (2017), Single-cell transcriptomes of human cell lines
- **Reference**: Li H, et al. *Nature Genetics* 2017;49:708-718
- **Objects**: `li$data` (expression matrix), `li$annotation` (cell annotations)
- **Obtaining**: Available from the original publication

#### `NCI60.RData`
- **Source**: NCI-60 Cancer Cell Line Panel
- **Reference**: Ross DT, et al. *Nature Genetics* 2000;24:227-235
- **Objects**: NCI60 dataset from the `ISLR` R package
- **Obtaining**: Available via `ISLR::NCI60`

#### `yan.RData`
- **Source**: Yan et al. (2013), Human preimplantation embryo scRNA-seq
- **Reference**: Yan L, et al. *Nature Structural & Molecular Biology* 2013;20:1131-1139
- **Accession**: GEO [GSE36552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36552)
- **Objects**: `yan$data` (expression matrix), `yan$annotation` (cell annotations)
- **Obtaining**: Available from GEO or the `scRNAseq` Bioconductor package

## Preparing Data

Example R code to prepare the datasets:

```r
# --- Khan dataset ---
# install.packages("ISLR")
library(ISLR)
data(Khan)
khan_mat   <- t(Khan$xtrain)  # genes x samples
khan_label <- Khan$ytrain
save(khan_mat, khan_label, file = "data/Khan.RData")

# --- NCI60 dataset ---
data(NCI60)
nci60_mat   <- t(NCI60$data)  # genes x samples
nci60_label <- NCI60$labs
save(nci60_mat, nci60_label, file = "data/NCI60.RData")

# --- Yan dataset ---
# BiocManager::install("scRNAseq")
# Or load from a preprocessed .RData file
yan_expr  <- t(yan$data)   # genes x cells
yan_label <- yan$annotation$type
save(yan_expr, yan_label, file = "data/yan.RData")

# --- Li dataset ---
# Load from a preprocessed .RData file
li_expr  <- li$data          # genes x cells
li_label <- li$annotation$type
save(li_expr, li_label, file = "data/li.RData")
```
