# IterativeV2

An improved version of our previous iterative clustering algorithm, with higher power DE testing and more flexibility.

```
devtools::install_github("mikejdeines/IterativeV2")
```

Performs Leiden clustering (Traag, V. A. et al. Sci Rep 9, 5233 (2019)) iteratively on Seurat objects to identify all clusters with a significant number of differentially-expressed genes.
The gene scoring metric was adapted from Tasic, B. et al. Nat Neurosci 19, 335–346 (2016), with minor modifications for differential expression testing and simplified cutoffs.

We are now using the weighted t-test DE functions from https://github.com/sergeyleikin/sc-sp-RNASeq. Please cite their paper if you use this package!

If a Seurat v5 object is being used, run JoinLayers() prior to clustering.

If the data has been normalized using SCTransform, run PrepSCTFindMarkers() prior to clustering.

This package was built and tested on Ubuntu 22.04.5 using Seurat v.5.3.0.
