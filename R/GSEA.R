#' GSEA as described in Subramanian et. al. 2005 PNAS
#' 
#' @param ges Gene expression signature (named numeric vector)
#' @param geneset Set of genes to test for enrichment in the GES (character vector)
#' @param score Power score for GSEA. Default value of 1 (standard for GSEA; use 0 for two-sample KS test)
#' @return An enrichment score for geneset in the GES.
GSEA <- function(ges,geneset,score = 1){
  
  # set the input values to local variables
  cur.geneset <- geneset
  cur.ges <- ges
  cur.score <- score
  
  # remove genes with missing values and zeros from the gene expression signature
  cur.ges <- cur.ges[which(!(cur.ges == 0 | is.na(cur.ges)))]
  
  # prune the gene set to include only those genes which are present in the gene expression signature
  cur.geneset <- cur.geneset[which(cur.geneset%in%names(cur.ges))]
  
  # sort the gene expression signature for most upregulated genes to most downregulated genes
  cur.sorted.ges <- sort(cur.ges,decreasing = TRUE)
  
  # match the genes in the gene set to the genes in the sorted gene expression signature
  cur.geneset.idx <- match((cur.geneset),names(cur.sorted.ges))
  
  # set running sum enrichment score values based on whether the gene is in the gene set or in the gene set's complement
  cur.geneset.enrichment <- rep((-1 / (length(cur.sorted.ges) - length(cur.geneset))),times = length(cur.sorted.ges))
  cur.geneset.enrichment[cur.geneset.idx] <- as.numeric(abs(cur.sorted.ges[cur.geneset.idx])^(cur.score)) / sum(as.numeric(abs(cur.sorted.ges[cur.geneset.idx])^(cur.score)))
  
  # compute the gene set running sum statistic
  cur.geneset.enrichment <- cumsum(cur.geneset.enrichment)
  
  # compute the gene set enrichment score as the supremum of the gene set running sum statistic
  cur.geneset.es <- cur.geneset.enrichment[which(abs(cur.geneset.enrichment) == max(abs(cur.geneset.enrichment)))]
  cur.geneset.es <- cur.geneset.es[1]
  
  # return the gene set enrichment score
  return(cur.geneset.es)
  
}