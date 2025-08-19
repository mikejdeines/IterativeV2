# IterativeV2

An improved version of our previous iterative clustering algorithm, with faster runtimes and more flexibility.\

```
devtools::install_github("mikejdeines/IterativeV2")
```

Performs Leiden clustering (Traag, V. A. et al. Sci Rep 9, 5233 (2019)) iteratively on Seurat objects to identify all clusters with a significant number of differentially-expressed genes.
The gene scoring metric was adapted from Tasic, B. et al. Nat Neurosci 19, 335–346 (2016), with minor modifications for differential expression testing and simplified cutoffs.\

This package requires the use of igraph's implementation of the Leiden algorithm.\
To use this implementation, create a conda environment using reticulate and install numpy 1.24.4, igraph, and leidenalg.
```
library(reticulate)
conda_create(envname = "scrnaseq")
use_condaenv("scrnaseq")
py_install("numpy==1.24.4")
py_install("igraph")
py_install("leidenalg")
```
If a Seurat v5 object is being used, run JoinLayers() prior to clustering.
