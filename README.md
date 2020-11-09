# NaRnEA

## NaRnEA (Non-parametric analytical Rank-based Enrichment Analysis)

Non-parametric analytical Rank-based Enrichment Analysis (`NaRnEA`) is a novel statistical method created by Griffin et al. in 2020 for the purpose of performing powerful and accurate gene set analysis. To read more about the theory behind NaRnEA and how it compares with Gene Set Enrichment Analysis (GSEA), which was created by Subramanian et al. in 2005, check out the preprint of the NaRnEA manuscript on bioRxiv ([link](https://www.biorxiv.org/search/NaRnEA)).

## Installing NaRnEA

NaRnEA is distributed free for scientific use as an R package; use the following code to download and install NaRnEA.
```{r}
devtools::install_github(repo = "califano-lab/NaRnEA", force = TRUE, build_vignettes = TRUE)
```

## Using NaRnEA
To learn more about how to use NaRnEA, check out the NaRnEA vignettes with the following code.
```{r}
library(NaRnEA)
browseVignettes(package = "NaRnEA")
```