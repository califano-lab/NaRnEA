#' Gene name convention map
#' 
#' A data frame containing gene names across five conventions: Ensembl transcripts (enst.values),
#' Ensembl genes (ensg.values), Entrez ID (entrez.values), Hugo gene names (hugo,values), and Biotype 
#' values (biotype.values). Contains 177,680 rows.
#' 
#' @docType data
#' 
#' @usage data(gene.name.map)
#' 
#' @format A data.frame
#' 
#' @source Generated using biomaRt
#' 
#' @examples 
#' data(gene.name.map)
#' \donttest{dim(gene.name.map)}
"gene.name.map"