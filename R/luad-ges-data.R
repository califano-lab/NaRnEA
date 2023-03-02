#' LUAD GES Vector
#' 
#' Gene expression signature calculated via tumor-vs-normal Mann-Whitney U-Test.
#' This object can be used to run the simplified vignette in the NaRnEA README.
#' For full instructions on how to generate this object, see the full LUAD vignette.
#' Genes in are labeled with Entrez IDs.
#'
#' @docType data
#' 
#' @usage data(LUAD_ges)
#' 
#' @format Two objects: `LUAD.regulon.list`, a list of regulon objects, each with parametrized network targets generated from ARACNe3;
#' and `LUAD.mwu.ges`, a vector of GES values for the tumor-vs-normal GES.
#' 
#' @source Objects generated from TCGA LUAD data.
#' 
#' @examples 
#' data(LUAD_ges)
#' \donttest{dim(LUAD_ges)}
"LUAD_ges"