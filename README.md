# pkgdown <img src="/man/Figures/Ghosn_Lab_logo_cropped.jpg" align="right"  height=200/>

# Data Matrix Normalization and Merging Strategies Minimize Batch-specific Systemic Variation in scRNA-Seq Data

This repository is intended to support reproducibility of our recent manuscript: "Data Matrix Normalization and Merging Strategies Minimize Batch-specific Systemic Variation in scRNA-Seq Data." The complete text and supplementary materials of the preprint are made freely available via bioRxiv, accessible [here](https://www.biorxiv.org/content/10.1101/2021.08.18.456898v1)

To enable rapid reproduction of the published figures, we have made a set of pre-assembled Seurat objects available as part of an R packaged "BatchNorm", available on github [here](https://github.com/Ghosn-Lab/BatchNorm).

![CMS_Overview](https://user-images.githubusercontent.com/50965273/136675082-4c9285c0-e1b2-45ab-a069-001e43a21c0c.jpg)

# Installation
To install the package (along with select published datasets in support of the manuscript) in `R`, run:

```r
if(!requireNamespace("devtools", quietly = TRUE)) {
 install.packages("devtools") 
}
devtools::install_github("Ghosn-Lab/BatchNorm")
```
The package should install within a few minutes, and all functions can be reproduced on a standard desktop or laptop without special hardware.

# Examples
To reproduce the results and figures as presented in the manuscript, you can follow along with our vignettes [here](https://ghosn-lab.github.io/BatchNorm/index.html)

## Biaxial Gating
[The first vignette](https://ghosn-lab.github.io/BatchNorm/articles/Biaxial_Gating.html) depicts the biaxial gating of a single sample (published as PBMC Sample 4-A). Following along with this vignette will reproduce the individual panels of Figure S3.



