#' Bulk RNA-Seq data from COAD
#' 
#' Bulk RNA-Seq data from The Cancer Genome Atlas (TCGA) for colon adenocarcinoma (COAD)
#' with HT-Seq counts data for 19,350 genes in 519 samples. Genes are in rows and samples 
#' are in columns. Genes are labeled with Entrez IDs.
#' 
#' @docType data
#' 
#' @usage data(TCGA_COAD)
#' 
#' @format A matrix.
#' 
#' @source Downloaded using TCGAbiolinks
#' 
#' @examples 
#' data(TCGA_COAD)
#' \donttest{dim(TCGA_COAD)}
"TCGA_COAD"