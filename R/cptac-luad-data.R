#' Protein Abundance data from LUAD
#' 
#' Protein abundance data downloaded from Clinical Proteomic Tumor Analysis Consortium (CPTAC)
#' for lung adenocarcinoma (LUAD). Data contains protein abundances for 10,316 proteins in 
#' 101 normal tissue samples and 110 tumor samples. Proteins are in rows and samples are in columns.
#' Protein are labeled with Hugo IDs.
#' 
#' @docType data
#' 
#' @usage data(CPTAC_LUAD)
#' 
#' @format A list with 'tissue' and 'tumor' elements, each a data.frame object.
#' 
#' @source Downloaded from CPTAC
#' 
#' @examples 
#' data(CPTAC_LUAD)
#' \donttest{dim(CPTAC_LUAD)}
"CPTAC_LUAD"