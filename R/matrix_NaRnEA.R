#' Runs the NaRnEA algorithm on multiple samples / regulons simultaneously.
#' 
#' @param ges.mat Matrix of gene expression signatures (genes X samples).
#' @param regulon.list List of regulon lists, with am and aw values.
#' @param seed.val Value of seed to ensure reproducibility. Default of 343.
#' @param min.targets Minimum number of targets needed between a regulon and a GES. Default of 30.
#' @return List of matrices; 'nes' and 'pes'.
#' @export
matrix_narnea <- function(ges.mat, regulon.list, seed.val = 343, min.targets = 30) {
  set.seed(seed.val)
  cat("Data prep...")
  
  ## remove edges not present in the NES; correct am values (nothing equal to 1, -1, or 0)
  regulon.list <- lapply(regulon.list, function(x) {
    # correct for targets
    ges.targets <- intersect(rownames(ges.mat), names(x$aw))
    aw.vec <- x$aw[ges.targets]
    am.vec <- x$am[ges.targets]
    # correct am values
    am.vec[which(am.vec == 1)] <- 0.999
    am.vec[which(am.vec == -1)] <- -0.999
    zero.samps <- which(am.vec == 0)
    am.vec[zero.samps] <- sample(c(0.001, -0.001), size = length(zero.samps), replace = TRUE)
    # reformat
    reg.list <- list('aw' = aw.vec, 'am' = am.vec)
    return(reg.list)
  })
  ## remove regulons with too few targets
  regulon.list <- regulon.list[which(sapply(regulon.list, function(x) {length(x$aw)}) >= min.targets)]
  if (length(regulon.list) == 0) {print('No regulons of adequate size'); return(NULL)}
  
  ## correct signature for zeros
  ges.mat <- apply(ges.mat, 2, function(x) {
    non.zero.num <- length(which(x == 0))
    if (non.zero.num == 0) { return(x) }
    pos.percent <- length(which(x > 0)) / non.zero.num
    neg.percent <- length(which(x < 0)) / non.zero.num
    non.zero.min <- min(abs(x[which(x != 0)]))
    x[which(x == 0)] <- runif(n = non.zero.num, min = 0, max = non.zero.min) *
      sample(c(-1, 1), size = non.zero.num, prob = c(neg.percent, pos.percent), replace = TRUE)
    return(x)
  })
  
  ## normalize GES
  R <- apply(ges.mat, 2, function(x) {rank(abs(x), ties.method = 'random')} )
  S <- apply(ges.mat, 2, function(x) {sign(x)} )
  ## get dimension parameters
  n <- ncol(ges.mat)
  g <- nrow(ges.mat)
  r <- length(regulon.list)
  ## create regulon matrices
  AM.mat <- matrix(0L, nrow = g, ncol = r)
  rownames(AM.mat) <- rownames(ges.mat); colnames(AM.mat) <- names(regulon.list)
  AW.mat <- matrix(0L, nrow = g, ncol = r)
  rownames(AW.mat) <- rownames(ges.mat); colnames(AW.mat) <- names(regulon.list)
  for (reg.name in names(regulon.list)) {
    AM.mat[names(regulon.list[[reg.name]]$am), reg.name] <- regulon.list[[reg.name]]$am
    AW.mat[names(regulon.list[[reg.name]]$aw), reg.name] <- regulon.list[[reg.name]]$aw
  }
  
  ## precalculated scalars
  E.r <- (g + 1) /2
  E.r2 <- (2*g^2 + 3*g + 1) / 6
  E.rs <- t(as.matrix((1 / g) * colSums(R * S)))
  ## precalculated matrices
  AM.abs.mat <- 1 - abs(AM.mat)
  AM.abs.mat[which(AM.mat == 0)] <- 0
  AW.AM.prod <- AW.mat * AM.mat
  AW.AM.abs.prod <- AW.mat * AM.abs.mat
  
  ## calculate NES
  cat("Calculating DES...")
  D.list <- directed_nes(R, S, AW.AM.prod, E.rs, E.r2, n)
  cat("Calculating UES...")
  U.list <- undirected_nes(R, S, AW.AM.abs.prod, E.r, E.r2, n)
  COV.nes <- nes_covariance(R, S, AW.AM.prod, AW.AM.abs.prod, E.r, E.rs, D.list$var, U.list$var, g)
  cat("Calculating NES...")
  NES.mat <- combine_nes(D.list$nes, U.list$nes, COV.nes)
  
  cat("Calculating PES...")
  ## calculate max D and U for each gene
  max.du <- lapply(names(regulon.list), function(x) {
    gene.order <- names(sort(regulon.list[[x]]$aw, decreasing = TRUE))
    aw.vec <- regulon.list[[x]]$aw[gene.order]
    abs.am.vec <- abs(regulon.list[[x]]$am[gene.order])
    ges.mag <- g:(g + 1 - length(gene.order))
    d.val <- sum(aw.vec * ges.mag * abs.am.vec)
    u.val <- sum(aw.vec * ges.mag * (1 - abs.am.vec))
    return(c(d.val, u.val))
  })
  max.du <- do.call('rbind', max.du)
  ## positive PES
  PES.pos.D <- matrix(rep(max.du[,1], n), nrow = r)
  PES.pos.U <- matrix(rep(max.du[,2], n), nrow = r)
  PES.pos.D.nes <- (PES.pos.D - D.list$exp) / sqrt(D.list$var)
  PES.pos.U.nes <- (PES.pos.U - U.list$exp) / sqrt(U.list$var)
  PES.pos.NES <- combine_nes(PES.pos.D.nes, PES.pos.U.nes, COV.nes)
  ## negative PES
  PES.neg.D <- PES.pos.D * (-1)
  PES.neg.U <- PES.pos.U
  PES.neg.D.nes <- (PES.neg.D - D.list$exp) / sqrt(D.list$var)
  PES.neg.U.nes <- (PES.neg.U - U.list$exp) / sqrt(U.list$var)
  PES.neg.NES <- combine_nes(PES.neg.D.nes, PES.neg.U.nes, COV.nes)
  ## combine, then normalize the NES scores
  pos.NES <- NES.mat > 0
  PES.comb.nes <- PES.pos.NES * pos.NES + PES.neg.NES * (!pos.NES)
  PES.mat <- NES.mat / abs(PES.comb.nes)
  
  cat("Done\n")
  return(list('NES' = NES.mat, 'PES' = PES.mat))
}

#' Computes the Directed NES
#' 
#' @param R Rank normalized ges matrix (genes x samples).
#' @param S Sign of ges matrix (genes x samples).
#' @param AW.AM.prod Product of the association weight (AW) and mode (AM) matrices created from the regulon list.
#' @param E.rs Expected value of r*s in each sample (1 x samples).
#' @param E.r2 Expected value of r^2 in each sample (scalar value).
#' @param n Number of samples (scalar value).
#' @return List of matrices, each (regulons x samples):
#' Enrichment score matrix 'es'; expected value 'exp'; variance 'var'; NES matrix 'NES'
#' @export
directed_nes <- function(R, S, AW.AM.prod, E.rs, E.r2, n) {
  D <- t(AW.AM.prod) %*% (R * S)
  ## calculate expected value
  D.e <- as.matrix(colSums(AW.AM.prod)) %*% E.rs
  ## calculate variance
  reg.squared.sum.product <- as.matrix(colSums((AW.AM.prod)**2))
  E.d.2 <- reg.squared.sum.product %*% (E.rs ** 2)
  E.d2 <- reg.squared.sum.product %*% matrix(rep(E.r2, n), nrow = 1)
  D.v <- E.d2 - E.d.2
  ## calculate NES
  D.nes <- (D - D.e) / sqrt(D.v)
  
  return(list('es' = D, 'exp' = D.e, 'var' = D.v, 'nes' = D.nes))
}

#' Computes the Undirected NES
#' 
#' @param R Rank normalized ges matrix (genes x samples).
#' @param S Sign of ges matrix (genes x samples).
#' @param AW.AM.abs.prod Product of the association weight (AW) and absolute value of the association mode (AM.abs) matrices created from the regulon list.
#' @param E.r Expected value of r in each sample (scalar value).
#' @param E.r2 Expected value of r^2 in each sample (scalar value).
#' @param n Number of samples (scalar value).
#' @return List of matrices, each (regulons x samples):
#' Enrichment score matrix 'es'; expected value 'exp'; variance 'var'; NES matrix 'NES'
#' @export
undirected_nes <- function(R, S, AW.AM.abs.prod, E.r, E.r2, n) {
  U <- t(AW.AM.abs.prod) %*% R
  ## calculate expected value
  U.e <- as.matrix(colSums(AW.AM.abs.prod)) %*% matrix(rep(E.r, n), nrow = 1)
  ## calculate variance
  reg.squared.sum.product.abs <- as.matrix(colSums((AW.AM.abs.prod)**2))
  E.u.2 <- reg.squared.sum.product.abs %*% matrix(rep(E.r**2, n), nrow = 1)
  E.u2 <- reg.squared.sum.product.abs %*% matrix(rep(E.r2, n), nrow = 1)
  U.v <- E.u2 - E.u.2
  ## calculate nes
  U.nes <- (U - U.e) / sqrt(U.v)
  
  return(list('es' = U, 'exp' = U.e, 'var' = U.v, 'nes' = U.nes))
}

#' Compute the Covariance of the Directed and Undirected NES.
#' 
#' @param R Rank normalized ges matrix (genes x samples).
#' @param S Sign of ges matrix (genes x samples).
#' @param AW.AM.prod Product of the association weight (AW) and mode (AM) matrices created from the regulon list.
#' @param AW.AM.abs.prod Product of the association weight (AW) and absolute value of the association mode (AM.abs) matrices created from the regulon list.
#' @param E.r Expected value of r in each sample (scalar value).
#' @param E.rs Expected value of r*s in each sample (1 x samples).
#' @param D.v Variance of Directed Enrichment Score (regulons x samples).
#' @param U.v Variance of Undirected Enrichment Score (regulons x samples).
#' @param g Number of genes in the gene expression signature matrix (scalar value).
#' @return Matrix of covariance values (regulons x samples).
#' @export
nes_covariance <- function(R, S, AW.AM.prod, AW.AM.abs.prod, E.r, E.rs, D.v, U.v, g) {
  cov.prod.mat <- as.matrix(colSums(AW.AM.prod * AW.AM.abs.prod))
  ## first component
  E.r2s <- as.matrix((1 / g) * colSums(R * R * S))
  E.du <- cov.prod.mat %*% t(E.r2s)
  ## second component
  E.d.u <- cov.prod.mat %*% E.rs * E.r
  ## combine
  COV.du <- E.du - E.d.u
  COV.nes <- COV.du / sqrt(D.v * U.v)
  
  return(COV.nes)
}

#' Generates a final NES from the combination of the Directed and Undirected NES.
#' 
#' @param D.nes NES scores for the directed enrichment (regulons x samples).
#' @param U.nes NES scores for the undirected enrichment (regulons x samples).
#' @param COV.nes Covariance of the undirected and directed enrichment NES (regulons x samples).
#' @return Matrix of integrated NES scores (regulons x samples).
#' @export
combine_nes <- function(D.nes, U.nes, COV.nes) {
  NES.pos <- (D.nes + U.nes) / sqrt(2 + 2*COV.nes)
  NES.neg <- (D.nes - U.nes) / sqrt(2 - 2*COV.nes)
  ## calculate p values
  p.pos <- pnorm(NES.pos, lower.tail = FALSE, log.p = TRUE)
  p.neg <- pnorm(NES.neg, lower.tail = TRUE, log.p = TRUE)
  ## combine p values
  p.dif <- (p.pos < p.neg)
  min.p <- (p.pos * p.dif + p.neg * (!p.dif))
  final.p <- min.p + log(2) + log1p(exp(min.p) / (-2))
  ## calculate final nes
  pos.nes <- qnorm(final.p - log(2), lower.tail = FALSE, log.p = TRUE)
  neg.nes <- qnorm(final.p - log(2), lower.tail = TRUE, log.p = TRUE)
  NES.mat <- (pos.nes * p.dif + neg.nes * (!p.dif))
  
  return(NES.mat)
}