# NaRnEA

## NaRnEA (Nonparametric analytical Rank-based Enrichment Analysis)

Nonparametric analytical Rank-based Enrichment Analysis (`NaRnEA`) is a novel statistical method created by Griffin et al. for the purpose of performing powerful and accurate gene set analysis. To read more about the theory behind `NaRnEA` and how it compares with Gene Set Enrichment Analysis (`GSEA`) and Virtual Inference of Protein-activity by Enriched Regulon analysis (`VIPER`), check out the preprint of the `NaRnEA` manuscript on bioRxiv ([link](https://www.biorxiv.org)).

## Installing the NaRnEA R package

To install the `NaRnEA` R package, download the zipped GitHub codebase, unzip it, and install locally with the `devtools` package using the code below (if the unzipped directory is named `NaRnEA-master`, change this to `NaRnEA` to ensure subsequent code functionality). The `NaRnEA` GitHub codebase also includes proteomic data and code for running `ARACNe3`, the newest implementation of the Algorithm for the Reconstruction of Accurate Cellular Networks, so that users can replicate the analysis from the `NaRnEA` manuscript. The code is freely available for academic research use; for more information, see the license associated with the `NaRnEA` GitHub codebase. To run `ARACNe3` users will need to install the latest version of JDK.

```{r}
devtools::install(pkg = "NaRnEA", build_vignettes = TRUE)
```

## Using NaRnEA and ARACNe3
To learn more about performing powerful and accurate gene set analysis with `NaRnEA`, check out the vignettes included in the `NaRnEA` GitHub codebase using the code below. These vignettes cover all analysis performed in the `NaRnEA` manuscript, including the reverse-engineering of transcriptional regulatory networks from gene expression data in The Cancer Genome Atlas (TCGA) for lung adenocarcinoma (LUAD), colon adenocarcinoma (COAD), and head-neck squamous cell carcinoma (HSNC) cancer subtypes with `ARACNe3`. 
To learn more about how to use `NaRnEA` for performing powerful and accurate gene set analysis, check out the `NaRnEA` vignettes with the following code.

```{r}
library(NaRnEA)
browseVignettes(package = "NaRnEA")
```
