#' Bulk RNA-Seq data from LUAD
#' 
#' Bulk RNA-Seq data from The Cancer Genome Atlas (TCGA) for lung adenocarcinoma (LUAD)
#' with HT-Seq counts data for 19,350 genes in 590 samples Genes are in rows and samples 
#' are in columns. Genes are labeled with Entrez IDs.
#' 
#' @docType data
#' 
#' @usage data(TCGA_LUAD)
#' 
#' @format A matrix.
#' 
#' @source Downloaded using TCGAbiolinks
#' 
#' @examples 
#' data(TCGA_LUAD)
#' \donttest(dim(TCGA_LUAD))
"TCGA_LUAD"