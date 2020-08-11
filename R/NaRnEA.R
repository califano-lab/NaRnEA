## NaRnEA = Non-parametric analytical Rank-based Enrichment Analysis
# ges = gene expression signature (named numeric vector)
# regulon = gene set parameterized with Regulation Confidence (likelihood) and Mode of Regulation (tfmode)
# minsize = minimum number of genes in the regulon that are also in the gene expression signature
# seed = random number generator seed to ensure that ties in the gene expression signature are broken randomly in a reproducible manner
# leading.edge = whether or not to compute the leading edge z-scores for the members of the gene set (logical) 

#' NaRnEA (Non-parameetric analytical rank-based enrichment analysis) function
#' 
#' @param ges Gene expression signature (named numeric vector)
#' @param regulon Gene set parameterized w/ Regulation Confidence (likelihood) and Mode of Regulation (tfmode)
#' @param minsize Minimum overlap between regulon and gene expression signature. Default of 30.
#' @param seed Randomly generated seed to ensure ties in the GES are broken in a reproducible manner. Default of 1.
#' @param leading.edge Flag to compute leading edge z-scores for the members of the gene set. False by default.
#' @return A list with the NES, PES, maximum and minimum possible NES, and the leading edge if specified.
NaRnEA <- function(ges, regulon, minsize = 30, seed = 1, leading.edge = FALSE){
  
  # set the seed
  set.seed(seed)
  
  # remove NAs and zeroes from the gene expression signature
  cur.ges <- ges[which(!(ges == 0 | is.na(ges)))]
  
  # subset the gene set to only include targets which are present in the gene expression signature
  cur.regul <- regulon
  cur.regul$tfmode <- regulon$tfmode[which(names(regulon$tfmode)%in%names(cur.ges))]
  cur.regul$likelihood <- regulon$likelihood[which(names(regulon$tfmode)%in%names(cur.ges))]
  
  # if the gene set is below the minimum size, return NA for all output values
  if(length(cur.regul$tfmode) < minsize){
    cur.res.list <- list(pes = NA, nes = NA, nes.max = NA, nes.min = NA, ledge = NA)
    return(cur.res.list)
  }
  
  # modify the Mode of Regulation values so that they are not equal to 1 or 0
  cur.regul$tfmode[which(abs(cur.regul$tfmode) == 1)] <- cur.regul$tfmode[which(abs(cur.regul$tfmode) == 1)]*(.999)
  cur.regul$tfmode[which(abs(cur.regul$tfmode) == 0)] <- (0.001)
  
  # compute the non-parametric gene expression signature transformation
  cur.ges.rank <- rank(abs(cur.ges), ties.method = "random")
  cur.ges.sign <- sign(cur.ges)
  
  # compute the number of genes in the gene expression signature
  cur.gene.num <- length(cur.ges)
  
  # identify the index of the gene set targets in the gene expression signature
  cur.idx <- match(names(cur.regul$tfmode),names(cur.ges))
  
  # compute the directed enrichment score for the gene set
  cur.directed.es.values <- (cur.regul$likelihood) * (cur.regul$tfmode) * (cur.ges.rank[cur.idx]) * (cur.ges.sign[cur.idx])
  cur.directed.es <- sum(cur.directed.es.values)
  
  # compute the expected value of the directed enrichment score for the gene set under the null hypothesis
  cur.directed.es.mean.values <- (cur.regul$likelihood) * (cur.regul$tfmode) * mean( (cur.ges.rank) * (cur.ges.sign))
  cur.directed.es.mean <- as.numeric(sum(cur.directed.es.mean.values))
  
  # compute the variance of the directed enrichment score for the gene set under the null hypothesis
  cur.directed.es.var.values <- (cur.regul$likelihood)^2 * (cur.regul$tfmode)^2 * (1/6) * (2*cur.gene.num^2 + 3*cur.gene.num + 1) - (cur.directed.es.mean.values)^2
  cur.directed.es.var <- as.numeric(sum(cur.directed.es.var.values))
  
  # compute the directed normalized enrichment score for the gene set
  cur.directed.nes <- (cur.directed.es - cur.directed.es.mean) / sqrt(cur.directed.es.var)
  
  # compute the undirected enrichment score for the gene set
  cur.undirected.es.values <- (cur.regul$likelihood) * (1 - abs(cur.regul$tfmode)) * (cur.ges.rank[cur.idx])
  cur.undirected.es <- as.numeric(sum(cur.undirected.es.values))
  
  # compute the expected value of the undirected enrichment score for the gene set under the null hypothesis
  cur.undirected.es.mean.values <- (cur.regul$likelihood) * (1 - abs(cur.regul$tfmode)) * (1/2) * (cur.gene.num + 1)
  cur.undirected.es.mean <- as.numeric(sum(cur.undirected.es.mean.values))
  
  # compute the variance of the undirected enrichment score for the gene set under the null hypothesis
  cur.undirected.es.var.values <- (cur.regul$likelihood)^2 * (1 - abs(cur.regul$tfmode))^2 * (1/6) * (2*cur.gene.num^2 + 3*cur.gene.num + 1) - (cur.undirected.es.mean.values)^2
  cur.undirected.es.var <- as.numeric(sum(cur.undirected.es.var.values))
  
  # compute the undirected normalized enrichment score for the gene set
  cur.undirected.nes <- (cur.undirected.es - cur.undirected.es.mean) / sqrt(cur.undirected.es.var)
  
  # compute the covariance between the directed and undirected enrichment scores for the gene set under the null hypothesis
  cur.directed.es.undirected.es.mean.values <- (cur.regul$likelihood)^2 * (cur.regul$tfmode) * (1 - abs(cur.regul$tfmode)) * mean( (cur.ges.rank)^2 * (cur.ges.sign))
  cur.directed.es.undirected.es.covar <- (cur.directed.es.mean.values) %*% t(cur.undirected.es.mean.values)
  diag(cur.directed.es.undirected.es.covar) <- cur.directed.es.undirected.es.mean.values
  cur.directed.es.undirected.es.covar <- sum(cur.directed.es.undirected.es.covar) - sum(cur.directed.es.mean.values) * sum(cur.undirected.es.mean.values)
  
  # compute the Pearson correlation between the directed and undirected enrichment scores for the gene set under the null hypothesis 
  cur.directed.es.undirected.es.cor <- (cur.directed.es.undirected.es.covar) / sqrt( (cur.directed.es.var) * (cur.undirected.es.var) )
  
  # compute the final normalized enrichment score for the gene set
  cur.nes.plus <- ( cur.directed.nes + cur.undirected.nes ) / sqrt( 2 + 2 * (cur.directed.es.undirected.es.cor) )
  cur.nes.minus <- ( cur.directed.nes - cur.undirected.nes ) / sqrt( 2 - 2 * (cur.directed.es.undirected.es.cor) )
  cur.plus.log.p.value <- pnorm(q = cur.nes.plus, lower.tail = FALSE, log.p = TRUE)
  cur.minus.log.p.value <- pnorm(q = cur.nes.minus, lower.tail = TRUE, log.p = TRUE)
  cur.min.log.p.value <- min(cur.plus.log.p.value, cur.minus.log.p.value)
  cur.final.log.p.value <- cur.min.log.p.value + log(2 - exp(cur.min.log.p.value))
  if(cur.plus.log.p.value == min(cur.plus.log.p.value, cur.minus.log.p.value)){
    cur.final.nes <- qnorm(p = (cur.final.log.p.value - log(2)), lower.tail = FALSE, log.p = TRUE)
  } else {
    cur.final.nes <- qnorm(p = (cur.final.log.p.value - log(2)), lower.tail = TRUE, log.p = TRUE)
  }
  
  # compute the maximum directed enrichment score for the gene set and corresponding undirected enrichment score for the gene set
  cur.directed.es.max.values <- sort( (cur.regul$likelihood) * (cur.regul$tfmode), decreasing = TRUE)
  cur.undirected.es.max.values <- (cur.regul$likelihood) * (1 - abs(cur.regul$tfmode))
  cur.undirected.es.max.values <- cur.undirected.es.max.values[match(names(cur.directed.es.max.values),names(cur.undirected.es.max.values))]
  if(sum(cur.directed.es.max.values > 0) == 0){
    cur.directed.es.max.ges.values <- as.numeric(rev(sort((cur.ges.rank*cur.ges.sign), decreasing = FALSE)[seq(from = 1, to = sum(cur.directed.es.max.values < 0), by = 1)]))
  } else if(sum(cur.directed.es.max.values < 0) == 0){
    cur.directed.es.max.ges.values <- as.numeric(sort((cur.ges.rank*cur.ges.sign), decreasing = TRUE)[seq(from = 1, to = sum(cur.directed.es.max.values > 0), by = 1)])
  } else {
    cur.directed.es.max.ges.values <- as.numeric(c(sort((cur.ges.rank*cur.ges.sign), decreasing = TRUE)[seq(from = 1, to = sum(cur.directed.es.max.values > 0), by = 1)],rev(sort((cur.ges.rank*cur.ges.sign), decreasing = FALSE)[seq(from = 1, to = sum(cur.directed.es.max.values < 0), by = 1)])))
  }
  cur.directed.es.max.values <- (cur.directed.es.max.values) * (cur.directed.es.max.ges.values)
  cur.undirected.es.max.values <- (cur.undirected.es.max.values) * abs(cur.directed.es.max.ges.values)
  cur.directed.es.max <- as.numeric(sum(cur.directed.es.max.values))
  cur.undirected.es.max <- as.numeric(sum(cur.undirected.es.max.values))
  
  # compute the minimum directed enrichment score for the gene set and corresponding undirected enrichment score for the gene set
  cur.directed.es.min.values <- sort( (cur.regul$likelihood) * (cur.regul$tfmode), decreasing = TRUE)
  cur.undirected.es.min.values <- (cur.regul$likelihood) * (1 - abs(cur.regul$tfmode))
  cur.undirected.es.min.values <- cur.undirected.es.min.values[match(names(cur.directed.es.min.values),names(cur.undirected.es.min.values))]
  if(sum(cur.directed.es.min.values > 0) == 0){
    cur.directed.es.min.ges.values <- as.numeric(rev(sort((cur.ges.rank*cur.ges.sign), decreasing = TRUE)[seq(from = 1, to = sum(cur.directed.es.min.values < 0), by = 1)]))
  } else if(sum(cur.directed.es.min.values < 0) == 0){
    cur.directed.es.min.ges.values <- as.numeric(sort((cur.ges.rank*cur.ges.sign), decreasing = FALSE)[seq(from = 1, to = sum(cur.directed.es.min.values > 0), by = 1)])
  } else {
    cur.directed.es.min.ges.values <- as.numeric(c(sort((cur.ges.rank*cur.ges.sign), decreasing = FALSE)[seq(from = 1, to = sum(cur.directed.es.min.values > 0), by = 1)],rev(sort((cur.ges.rank*cur.ges.sign), decreasing = TRUE)[seq(from = 1, to = sum(cur.directed.es.min.values < 0), by = 1)])))
  }
  cur.directed.es.min.values <- (cur.directed.es.min.values) * (cur.directed.es.min.ges.values)
  cur.undirected.es.min.values <- (cur.undirected.es.min.values) * abs(cur.directed.es.min.ges.values)
  cur.directed.es.min <- as.numeric(sum(cur.directed.es.min.values))
  cur.undirected.es.min <- as.numeric(sum(cur.undirected.es.min.values))
  
  # compute the maximum directed normalized enrichment score for the gene set
  cur.directed.nes.max <- (cur.directed.es.max - cur.directed.es.mean) / sqrt(cur.directed.es.var)
  
  # compute the maximum undirected normalized enrichment score for the gene set
  cur.undirected.nes.max <- (cur.undirected.es.max - cur.undirected.es.mean) / sqrt(cur.undirected.es.var)
  
  # compute the minimum directed normalized enrichment score for the gene set
  cur.directed.nes.min <- (cur.directed.es.min - cur.directed.es.mean) / sqrt(cur.directed.es.var)
  
  # compute the minimum undirected normalized enrichment score for the gene set
  cur.undirected.nes.min <- (cur.undirected.es.min - cur.undirected.es.mean) / sqrt(cur.undirected.es.var)
  
  # compute the maximum final normalized enrichment score for the gene set
  cur.nes.plus <- ( cur.directed.nes.max + cur.undirected.nes.max ) / sqrt( 2 + 2 * (cur.directed.es.undirected.es.cor) )
  cur.nes.minus <- ( cur.directed.nes.max - cur.undirected.nes.max ) / sqrt( 2 - 2 * (cur.directed.es.undirected.es.cor) )
  cur.plus.log.p.value <- pnorm(q = cur.nes.plus, lower.tail = FALSE, log.p = TRUE)
  cur.minus.log.p.value <- pnorm(q = cur.nes.minus, lower.tail = TRUE, log.p = TRUE)
  cur.min.log.p.value <- min(cur.plus.log.p.value, cur.minus.log.p.value)
  cur.final.log.p.value <- cur.min.log.p.value + log(2 - exp(cur.min.log.p.value))
  cur.final.nes.max <- qnorm(p = (cur.final.log.p.value - log(2)), lower.tail = FALSE, log.p = TRUE)
  
  # compute the minimum final normalized enrichment score for the gene set
  cur.nes.plus <- ( cur.directed.nes.min + cur.undirected.nes.min ) / sqrt( 2 + 2 * (cur.directed.es.undirected.es.cor) )
  cur.nes.minus <- ( cur.directed.nes.min - cur.undirected.nes.min ) / sqrt( 2 - 2 * (cur.directed.es.undirected.es.cor) )
  cur.plus.log.p.value <- pnorm(q = cur.nes.plus, lower.tail = FALSE, log.p = TRUE)
  cur.minus.log.p.value <- pnorm(q = cur.nes.minus, lower.tail = TRUE, log.p = TRUE)
  cur.min.log.p.value <- min(cur.plus.log.p.value, cur.minus.log.p.value)
  cur.final.log.p.value <- cur.min.log.p.value + log(2 - exp(cur.min.log.p.value))
  cur.final.nes.min <- (-1) * qnorm(p = (cur.final.log.p.value - log(2)), lower.tail = FALSE, log.p = TRUE)
  
  # compute the final proportional enrichment score for the gene set
  if(is.na(cur.final.nes)){
    cur.res.list <- list(nes = NA, pes = NA)
    return(cur.res.list)
  } else if(cur.final.nes > 0){
    cur.final.pes <- (cur.final.nes) / abs(cur.final.nes.max)
  } else {
    cur.final.pes <- (cur.final.nes) / abs(cur.final.nes.min)
  }
  
  # computing the leading edge for the gene set if desired
  if(leading.edge){
    if(cur.final.pes > 0){
      cur.ledge.values <- cur.directed.es.values + cur.undirected.es.values
      cur.ledge.p.values <- sapply(names(cur.ledge.values),function(cur.gene){
        cur.gene.ic.value <- as.numeric(cur.regul$likelihood[match(cur.gene,names(cur.regul$tfmode))])
        cur.gene.mor.value <- as.numeric(cur.regul$tfmode[match(cur.gene,names(cur.regul$tfmode))])
        cur.gene.null.ledge.values <- (cur.gene.ic.value)*(cur.gene.mor.value)*(cur.ges.rank)*(cur.ges.sign) + (cur.gene.ic.value)*(1 - abs(cur.gene.mor.value))*(cur.ges.rank)
        cur.gene.ledge.p.value <- (1 + sum((cur.gene.null.ledge.values > as.numeric(cur.ledge.values[match(cur.gene,names(cur.ledge.values))]))))/(1 + length(cur.gene.null.ledge.values))
        return(cur.gene.ledge.p.value)
      })
      cur.ledge.z.scores <- qnorm(p = cur.ledge.p.values, lower.tail = FALSE)
    } else {
      cur.ledge.values <- cur.directed.es.values - cur.undirected.es.values
      cur.ledge.p.values <- sapply(names(cur.ledge.values),function(cur.gene){
        cur.gene.ic.value <- as.numeric(cur.regul$likelihood[match(cur.gene,names(cur.regul$tfmode))])
        cur.gene.mor.value <- as.numeric(cur.regul$tfmode[match(cur.gene,names(cur.regul$tfmode))])
        cur.gene.null.ledge.values <- (cur.gene.ic.value)*(cur.gene.mor.value)*(cur.ges.rank)*(cur.ges.sign) - (cur.gene.ic.value)*(1 - abs(cur.gene.mor.value))*(cur.ges.rank)
        cur.gene.ledge.p.value <- (1 + sum((cur.gene.null.ledge.values < as.numeric(cur.ledge.values[match(cur.gene,names(cur.ledge.values))]))))/(1 + length(cur.gene.null.ledge.values))
        return(cur.gene.ledge.p.value)
      })
      cur.ledge.z.scores <- qnorm(p = cur.ledge.p.values, lower.tail = TRUE)
    }
    cur.final.ledge <- cur.ledge.z.scores
  } else {
    cur.final.ledge <- NA
  }
  
  # return the normalized enrichment score, the proportional enrichment score, the maximum normalized enrichment score, the minimum normalized enrichment score, and the leading edge (if desired) for the gene set
  cur.res.list <- list(pes = cur.final.pes, nes = cur.final.nes, nes.max = cur.final.nes.max, nes.min = cur.final.nes.min, ledge = cur.final.ledge)
  return(cur.res.list)
  
}
