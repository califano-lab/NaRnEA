#' Bulk RNA-Seq data from HNSC
#' 
#' Bulk RNA-Seq data from The Cancer Genome Atlas (TCGA) for head and neck squamous cell carcinoma (HNSC)
#' with HT-Seq counts data for 19,350 genes in 543 samples. Genes are in rows and samples 
#' are in columns. Genes are labeled with Entrez IDs.
#' 
#' @docType data
#' 
#' @usage data(TCGA_HNSC)
#' 
#' @format A matrix.
#' 
#' @source Downloaded using TCGAbiolinks
#' 
#' @examples 
#' data(TCGA_HNSC)
#' \donttest{dim(TCGA_HNSC)}
"TCGA_HNSC"