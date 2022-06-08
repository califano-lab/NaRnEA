#' Protein Abundance data from HNSC
#' 
#' Protein abundance data downloaded from Clinical Proteomic Tumor Analysis Consortium (CPTAC)
#' for head and neck squamous cell carcinoma (HNSC). Data contains protein abundances for 9,666 proteins in 
#' 63 normal tissue samples and 109 tumor samples. Proteins are in rows and samples are in columns.
#' Protein are labeled with Hugo IDs.
#' 
#' @docType data
#' 
#' @usage data(CPTAC_HNSC)
#' 
#' @format A list with 'tissue' and 'tumor' elements, each a data.frame object.
#' 
#' @source Downloaded from CPTAC
#' 
#' @examples 
#' data(CPTAC_HNSC)
#' \donttest{dim(CPTAC_HNSC)}
"CPTAC_HNSC"