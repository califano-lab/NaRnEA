#' Protein Abundance data from COAD
#' 
#' Protein abundance data downloaded from Clinical Proteomic Tumor Analysis Consortium (CPTAC)
#' for colon adenocarcinoma (COAD). Data contains protein abundances for 8,067 proteins in 
#' 100 normal tissue samples and 97 tumor samples. Proteins are in rows and samples are in columns.
#' Protein are labeled with Hugo IDs.
#' 
#' @docType data
#' 
#' @usage data(CPTAC_COAD)
#' 
#' @format A list with 'tissue' and 'tumor' elements, each a data.frame object.
#' 
#' @source Downloaded from CPTAC
#' 
#' @examples 
#' data(CPTAC_COAD)
#' \donttest(dim(CPTAC_COAD))
"CPTAC_COAD"