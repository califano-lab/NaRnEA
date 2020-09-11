# NaRnEA

## NaRnEA (Non-parametric analytical rank-based enrichment analysis)

Non-parametric analytical Rank-based Enrichment Analysis (NaRnEA) enables users to perform accurate and powerful gene set analysis in a non-parametric, frequentist manner.

## Installing NaRnEA

NaRnEA is distributed free for scientific use as an R package; use the following code to download and install NaRnEA.
```{r}
devtools::install_github(repo = "califano-lab/NaRnEA", force = TRUE, build_vignettes = TRUE)
```

## NaRnEA Theory and Walkthroughs
To learn more about using NaRnEA and the theory behind NaRnEA, check out the NaRnEA vignettes.
```{r}
library(NaRnEA)
browseVignettes(package = "NaRnEA")
```