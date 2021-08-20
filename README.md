# NaRnEA

## NaRnEA (Non-parametric analytical Rank-based Enrichment Analysis)

Non-parametric analytical Rank-based Enrichment Analysis (`NaRnEA`) is a novel statistical method created by Griffin et al. for the purpose of performing powerful and accurate gene set analysis. To read more about the theory behind `NaRnEA` and how it compares with Gene Set Enrichment Analysis (GSEA) and Virtual Inference of Protein-activity by Enriched Regulon analysis (VIPER), check out the preprint of the `NaRnEA` manuscript on bioRxiv ([link](https://www.biorxiv.org)).

## Installing the NaRnEA R package

The `NaRnEA` R package can be downloaded with the following code and is distributed free for scientific use.
```{r}
devtools::install_github(repo = "califano-lab/NaRnEA", force = TRUE, build_vignettes = TRUE)
```

The NaRnEA R package also includes code necessary for running the Algorithm for the Reconstruction of Accurate Cellular Networks version 3 (`ARACNe3`).  

## Using NaRnEA and ARACNe3
To learn more about how to use `NaRnEA` for performing powerful and accurate gene set analysis in a non-parametric manner, check out the `NaRnEA` vignettes with the following code.
```{r}
library(NaRnEA)
browseVignettes(package = "NaRnEA")
```

