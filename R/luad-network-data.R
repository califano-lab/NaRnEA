#' LUAD Network List
#' 
#' Network list generated from tumor samples using ARACNe3.
#' This object can be used to run the simplified vignette in the NaRnEA README.
#' For full instructions on how to generate this object, see the full LUAD vignette.
#' Genes are labeled with Entrez IDs.
#'
#' @docType data
#' 
#' @usage data(LUAD_network)
#' 
#' @format Two objects: `LUAD.regulon.list`, a list of regulon objects, each with parametrized network targets generated from ARACNe3;
#' and `LUAD.mwu.ges`, a vector of GES values for the tumor-vs-normal GES.
#' 
#' @source Objects generated from TCGA LUAD data.
#' 
#' @examples 
#' data(LUAD_network)
#' \donttest{dim(LUAD_network)}
"LUAD_network"