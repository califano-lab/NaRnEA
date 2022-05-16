#' Nonparametric analytical rank-based enrichment analysis (NaRnEA) function
#' 
#' @param signature Gene expression signature (named numeric vector)
#' @param association.weight Association Weight values for gene set members (named numeric vector)
#' @param association.mode Association Mode values for gene set members (named numeric vector)
#' @param ledge Flag to compute ledge p-values for gene set members.
#' @param minimum.size Minimum number of gene set members. Default of 30.
#' @param seed Random number generator seed to ensure reproducibility manner. Default of 1.
#' @return A list with the Proportional Enrichment Score (PES), Normalized Enrichment Score (NES), and the leading edge p-values (NA if ledge set to FALSE).
#' @export

NaRnEA <- function(signature, association.weight, association.mode, ledge = TRUE, minimum.size = 30, seed = 1){
  
  # set the seed
  set.seed(seed)
  
  # set the input values to local variables
  cur.sig <- signature
  cur.rc.values <- association.weight
  cur.mor.values <- association.mode
  
  # check that the Regulation Confidence values are positive
  if(!prod(cur.rc.values > 0)){
    stop("... regulation.confidence values must be positive ...")
  } 
  
  # check that the Mode of Regulation values are between (-1) and (1)
  if(!prod(abs(cur.mor.values) <= 1)){
    stop("... mode.of.regulation values must be between (-1) and (1) ...")
  }
  
  # correct and Mode of Regulation values which are exactly equal to 0 or exactly equal to 1 
  cur.mor.values[which(abs(cur.mor.values) == 1)] <- 0.999*sign(cur.mor.values[which(abs(cur.mor.values) == 1)])
  cur.mor.values[which(cur.mor.values == 0)] <- sample(c(0.001, -0.001), size = length(which(cur.mor.values == 0)), replace = TRUE, prob = c(mean(cur.mor.values[which(cur.mor.values != 0)] > 0), mean(cur.mor.values[which(cur.mor.values != 0)] < 0)))
  
  # check that the names for the Regulation Confidence and Mode of Regulation match and return an error message if they do not
  if(identical(names(cur.rc.values),names(cur.mor.values))){
    cur.target.values <- names(cur.rc.values)
  } else {
    stop("... regulation.confidence and mode.of.regulation must have the same names ...")
  }
  
  # compute the nonparametric signature transformation if no missing values or repeated entries are present; if missing values or repeated entries are present in the signature, return an error message
  if(anyNA(cur.sig)){
    stop("... signature cannot contain missing values ...")
  } else if(!prod(!duplicated(names(cur.sig)))){
    stop("... signature cannot contain duplicate entries ...")
  } else if(prod(cur.sig != 0)){
    cur.np.sig <- rank(x = abs(as.numeric(cur.sig)), ties.method = "random")*sign(as.numeric(cur.sig))
    names(cur.np.sig) <- names(cur.sig)
  } else {
  	set.seed(seed)
	cur.sig[which(cur.sig == 0)] <- runif(n = sum(cur.sig == 0), min = 0, max = min(abs(cur.sig[which(cur.sig != 0)])))*sample(c(-1,1), size = sum(cur.sig == 0), replace = TRUE, prob = c(mean(cur.sig[which(cur.sig != 0)] < 0),mean(cur.sig[which(cur.sig != 0)] > 0)))
    cur.np.sig <- rank(x = abs(as.numeric(cur.sig)), ties.method = "random")*sign(as.numeric(cur.sig))
    names(cur.np.sig) <- names(cur.sig)
  }
  
  # check that all of the target names are in the gene expression signature and return an error message if any are not
  if(mean(cur.target.values%in%names(cur.np.sig)) < 1){
    stop("... not all of the gene set members are present in the gene expression signature ...")
  }
  
  # return an error message if the gene set is below the minimum size; if not, compute the indices for each gene set member
  if(length(cur.target.values) < minimum.size){
    stop("... gene set is too small ...")
  } else {
    cur.idx.values <- match(cur.target.values,names(cur.np.sig))
  }
  
  # compute the Directed Enrichment Score for the gene set
  cur.directed.es.values <- (cur.rc.values)*(cur.mor.values)*(cur.np.sig[cur.idx.values]) 
  cur.directed.es <- sum(cur.directed.es.values)
  
  # compute the Normalized Directed Enrichment Score for the gene set
  cur.directed.es.mean.values <- (cur.rc.values)*(cur.mor.values)*mean(cur.np.sig)
  cur.directed.es.mean <- sum(cur.directed.es.mean.values)
  
  cur.directed.es.var.values <- (cur.rc.values^2)*(cur.mor.values^2)*(2*length(cur.np.sig)^2 + 3*length(cur.np.sig) + 1)*(1/6) - (cur.directed.es.mean.values^2)
  cur.directed.es.var <- sum(cur.directed.es.var.values)
  
  cur.directed.nes <- (cur.directed.es - cur.directed.es.mean)/sqrt(cur.directed.es.var)
  
  # compute the Undirected Enrichment Score for the gene set
  cur.undirected.es.values <- (cur.rc.values)*(1 - abs(cur.mor.values))*(abs(cur.np.sig)[cur.idx.values])
  cur.undirected.es <- sum(cur.undirected.es.values)
  
  # compute the Normalized Undirected Enrichment Score for the gene set
  cur.undirected.es.mean.values <- (cur.rc.values)*(1 - abs(cur.mor.values))*(length(cur.np.sig) + 1)*(1/2)
  cur.undirected.es.mean <- sum(cur.undirected.es.mean.values)
  
  cur.undirected.es.var.values <- (cur.rc.values^2)*((1 - abs(cur.mor.values))^2)*(2*length(cur.np.sig)^2 + 3*length(cur.np.sig) + 1)*(1/6) - (cur.undirected.es.mean.values^2)
  cur.undirected.es.var <- sum(cur.undirected.es.var.values)
  
  cur.undirected.nes <- (cur.undirected.es - cur.undirected.es.mean)/sqrt(cur.undirected.es.var)
  
  # compute the covariance between the Normalized Directed Enrichment Score and the Normalized Undirected Enrichment Score under the null hypothesis
  cur.directed.es.undirected.es.mean.mat <- matrix(data = cur.directed.es.mean.values, nrow = length(cur.directed.es.mean.values), ncol = length(cur.directed.es.mean.values), byrow = TRUE)*matrix(data = cur.undirected.es.mean.values, nrow = length(cur.undirected.es.mean.values), ncol = length(cur.undirected.es.mean.values), byrow = FALSE)
  
  diag(cur.directed.es.undirected.es.mean.mat) <- (cur.rc.values^2)*(cur.mor.values)*(1 - abs(cur.mor.values))*mean(cur.np.sig*abs(cur.np.sig))
  
  cur.directed.es.undirected.es.mean <- sum(cur.directed.es.undirected.es.mean.mat)
  
  cur.directed.es.undirected.es.covariance <- cur.directed.es.undirected.es.mean - (cur.directed.es.mean)*(cur.undirected.es.mean)
  
  cur.directed.nes.undirected.nes.covariance <- cur.directed.es.undirected.es.covariance/sqrt((cur.directed.es.var)*(cur.undirected.es.var))
  
  # compute the final Normalized Enrichment Score for the gene set
  cur.nes.plus <- (cur.directed.nes + cur.undirected.nes)/sqrt(2 + 2*(cur.directed.nes.undirected.nes.covariance))
  cur.plus.log.p.value <- pnorm(q = cur.nes.plus, lower.tail = FALSE, log.p = TRUE)
  
  cur.nes.minus <- (cur.directed.nes - cur.undirected.nes)/sqrt(2 - 2*(cur.directed.nes.undirected.nes.covariance))
  cur.minus.log.p.value <- pnorm(q = cur.nes.minus, lower.tail = TRUE, log.p = TRUE)
  
  cur.final.log.p.value <- (min(c(cur.plus.log.p.value,cur.minus.log.p.value))) + log(2) + log1p( (exp( min(c(cur.plus.log.p.value,cur.minus.log.p.value)) )/(-2) ) )
  
  if(cur.plus.log.p.value < cur.minus.log.p.value){
    cur.final.nes <- qnorm(p = (cur.final.log.p.value - log(2)), lower.tail = FALSE, log.p = TRUE)
  } else {
    cur.final.nes <- qnorm(p = (cur.final.log.p.value - log(2)), lower.tail = TRUE, log.p = TRUE)
  }
  
  # compute the Proportional Enrichment Score for the gene set
  if(cur.final.nes > 0){
  
    # compute the target indices which correspond with supreme positive gene set enrichment
    max.idx.values <- cur.idx.values
    if(sum(cur.mor.values > 0) > 0){
      # max.idx.values[which(cur.mor.values > 0)] <- match(names(rev(sort(cur.np.sig[which(cur.np.sig > 0)], decreasing = TRUE)[seq(from = 1, to = length(which(cur.mor.values > 0)), by = 1)])[rank(cur.rc.values[which(cur.mor.values > 0)], ties.method = "random")]), names(cur.np.sig))
      max.idx.values[which(cur.mor.values > 0)] <- match(names(rev(sort(cur.np.sig, decreasing = TRUE)[seq(from = 1, to = length(which(cur.mor.values > 0)), by = 1)])[rank(cur.rc.values[which(cur.mor.values > 0)], ties.method = "random")]), names(cur.np.sig))
    }
    if(sum(cur.mor.values < 0) > 0){
      # max.idx.values[which(cur.mor.values < 0)] <- match(names(rev(sort(cur.np.sig[which(cur.np.sig < 0)], decreasing = FALSE)[seq(from = 1, to = length(which(cur.mor.values < 0)), by = 1)])[rank(cur.rc.values[which(cur.mor.values < 0)], ties.method = "random")]), names(cur.np.sig))
      max.idx.values[which(cur.mor.values < 0)] <- match(names(rev(sort(cur.np.sig, decreasing = FALSE)[seq(from = 1, to = length(which(cur.mor.values < 0)), by = 1)])[rank(cur.rc.values[which(cur.mor.values < 0)], ties.method = "random")]), names(cur.np.sig))
    }
    
    # compute the Directed Enrichment Score for the supreme positive gene set enrichment target indices
    max.directed.es.values <- (cur.rc.values)*(cur.mor.values)*(cur.np.sig[max.idx.values])
    max.directed.es <- sum(max.directed.es.values)
    
    # compute the Normalized Directed Enrichment Score for the supreme positive gene set enrichment target indices
    max.directed.nes <- (max.directed.es - cur.directed.es.mean)/sqrt(cur.directed.es.var)
    
    # compute the Undirected Enrichment Score for the supreme positive gene set enrichment target indices
    max.undirected.es.values <- (cur.rc.values)*(1 - abs(cur.mor.values))*(abs(cur.np.sig)[max.idx.values])
    max.undirected.es <- sum(max.undirected.es.values)
    
    # compute the Normalized Undirected Enrichment Score for the supreme positive gene set enrichment target indices
    max.undirected.nes <- (max.undirected.es - cur.undirected.es.mean)/sqrt(cur.undirected.es.var)
    
    # compute the final Normalized Enrichment Score for the supreme positive gene set enrichment target indices
    max.nes.plus <- (max.directed.nes + max.undirected.nes)/sqrt(2 + 2*(cur.directed.nes.undirected.nes.covariance))
    max.plus.log.p.value <- pnorm(q = max.nes.plus, lower.tail = FALSE, log.p = TRUE)
    
    max.final.log.p.value <- (max.plus.log.p.value) + log(2) + log1p((exp(max.plus.log.p.value)/(-2)))
    
    max.final.nes <- qnorm(p = (max.final.log.p.value - log(2)), lower.tail = FALSE, log.p = TRUE)
    
    # compute the Proportional Enrichment Score for the gene set
    cur.final.pes <- (cur.final.nes/abs(max.final.nes))
    cur.final.pes <- min(c(cur.final.pes,1))
    
  } else {
  
    # compute the target indices which correspond with supreme negative gene set enrichment
    min.idx.values <- cur.idx.values
    if(sum(cur.mor.values > 0) > 0){
      # min.idx.values[which(cur.mor.values > 0)] <- match(names(rev(sort((-1*cur.np.sig)[which((-1*cur.np.sig) > 0)], decreasing = TRUE)[seq(from = 1, to = length(which(cur.mor.values > 0)), by = 1)])[rank(cur.rc.values[which(cur.mor.values > 0)], ties.method = "random")]), names((-1*cur.np.sig)))
      min.idx.values[which(cur.mor.values > 0)] <- match(names(rev(sort((-1*cur.np.sig), decreasing = TRUE)[seq(from = 1, to = length(which(cur.mor.values > 0)), by = 1)])[rank(cur.rc.values[which(cur.mor.values > 0)], ties.method = "random")]), names((-1*cur.np.sig)))
    }
    if(sum(cur.mor.values < 0) > 0){
      # min.idx.values[which(cur.mor.values < 0)] <- match(names(rev(sort((-1*cur.np.sig)[which((-1*cur.np.sig) < 0)], decreasing = FALSE)[seq(from = 1, to = length(which(cur.mor.values < 0)), by = 1)])[rank(cur.rc.values[which(cur.mor.values < 0)], ties.method = "random")]), names((-1*cur.np.sig)))
      min.idx.values[which(cur.mor.values < 0)] <- match(names(rev(sort((-1*cur.np.sig), decreasing = FALSE)[seq(from = 1, to = length(which(cur.mor.values < 0)), by = 1)])[rank(cur.rc.values[which(cur.mor.values < 0)], ties.method = "random")]), names((-1*cur.np.sig)))
    }
    
    # compute the Directed Enrichment Score for the supreme negative gene set enrichment target indices
    min.directed.es.values <- (cur.rc.values)*(cur.mor.values)*(cur.np.sig[min.idx.values])
    min.directed.es <- sum(min.directed.es.values)
    
    # compute the Normalized Directed Enrichment Score for the supreme negative gene set enrichment target indices
    min.directed.nes <- (min.directed.es - cur.directed.es.mean)/sqrt(cur.directed.es.var)
    
    # compute the Undirected Enrichment Score for the supreme negative gene set enrichment target indices
    min.undirected.es.values <- (cur.rc.values)*(1 - abs(cur.mor.values))*(abs(cur.np.sig)[min.idx.values])
    min.undirected.es <- sum(min.undirected.es.values)
    
    # compute the Normalized Undirected Enrichment Score for the supreme negative gene set enrichment target indices
    min.undirected.nes <- (min.undirected.es - cur.undirected.es.mean)/sqrt(cur.undirected.es.var)
    
    # compute the final Normalized Enrichment Score for the supreme negative gene set enrichment target indices
    min.nes.minus <- (min.directed.nes - min.undirected.nes)/sqrt(2 - 2*(cur.directed.nes.undirected.nes.covariance))
    min.minus.log.p.value <- pnorm(q = min.nes.minus, lower.tail = TRUE, log.p = TRUE)
    
    min.final.log.p.value <- (min.minus.log.p.value) + log(2) + log1p((exp(min.minus.log.p.value)/(-2)))
    
    min.final.nes <- qnorm(p = (min.final.log.p.value - log(2)), lower.tail = TRUE, log.p = TRUE)
    
    # compute the Proportional Enrichment Score for the gene set
    cur.final.pes <- (cur.final.nes/abs(min.final.nes))
    cur.final.pes <- max(c(cur.final.pes,-1))
  
  }
  
  # compute the ledge p-values for the gene set if desired
  if(ledge){
    if(cur.final.pes > 0){
      
      cur.ledge.p.values <- sapply(cur.target.values,function(sub.target.value){
        
        # compute the ledge score for the current gene set member
        sub.ledge.score <- (1 - abs(cur.mor.values))[match(sub.target.value, names(cur.mor.values))]*(abs(cur.np.sig)[cur.idx.values[match(sub.target.value, names(cur.mor.values))]]) + (cur.mor.values)[match(sub.target.value, names(cur.mor.values))]*(cur.np.sig[cur.idx.values[match(sub.target.value, names(cur.mor.values))]])
        
        # compute the null ledge scores for the current gene set member
        sub.ledge.null <- (1 - abs(cur.mor.values))[match(sub.target.value, names(cur.mor.values))]*(abs(cur.np.sig)) + (cur.mor.values)[match(sub.target.value, names(cur.mor.values))]*(cur.np.sig)
        
        # compute the ledge p-value for the current gene set member
        sub.ledge.p.value <- mean((sub.ledge.null >= sub.ledge.score))
        return(sub.ledge.p.value)
        
      })
      
    } else {
      
      cur.ledge.p.values <- sapply(cur.target.values,function(sub.target.value){
        
        # compute the ledge score for the current gene set member
        sub.ledge.score <- (1 - abs(cur.mor.values))[match(sub.target.value, names(cur.mor.values))]*(abs(cur.np.sig)[cur.idx.values[match(sub.target.value, names(cur.mor.values))]]) - (cur.mor.values)[match(sub.target.value, names(cur.mor.values))]*(cur.np.sig[cur.idx.values[match(sub.target.value, names(cur.mor.values))]])
        
        # compute the null ledge scores for the current gene set member
        sub.ledge.null <- (1 - abs(cur.mor.values))[match(sub.target.value, names(cur.mor.values))]*(abs(cur.np.sig)) - (cur.mor.values)[match(sub.target.value, names(cur.mor.values))]*(cur.np.sig)
        
        # compute the ledge p-value for the current gene set member
        sub.ledge.p.value <- mean((sub.ledge.null >= sub.ledge.score))
        return(sub.ledge.p.value)
        
      })
      
    }
  } else {
    cur.ledge.p.values <- rep(NA, times = length(cur.target.values))
    names(cur.ledge.p.values) <- cur.target.values
  }
  
  # return the Proportional Enrichment Score, the final Normalized Enrichment Score, and the ledge p-values (if desired) for the gene set 
  cur.final.res.list <- list(pes = cur.final.pes, nes = cur.final.nes, ledge = cur.ledge.p.values)
  return(cur.final.res.list)

}
