#' Nonparametric analytical Rank-based Enrichment Analysis plotting function
#' 
#' @param signature Gene expression signature (named numeric vector)
#' @param association.weight Association Weight values for gene set members (named numeric vector)
#' @param association.mode Association Mode values for gene set members (named numeric vector)
#' @param minimum.size Minimum number of gene set members. Default of 30.
#' @param seed Random number generator seed to ensure reproducibility manner. Default of 1.
#' @param title.text Text for title of the plot. 
#' @param x.label.text Text for x-axis label of the plot. 
#' @param line.width Line width for the gene set member segments in the nonparametric gene expression signature. 
#' @param lower.color Color for the gene set members with negative Association Mode values.
#' @param middle.color Color for the gene set members with Association Mode values close to zero.
#' @param upper.color Color for the gene set members with positive Association Mode values. 
#' @param conf.int Flag to compute the confidence interval for the NaRnEA Proportional Enrichment Score using a Fisher Transformation and gene set member sampling with replacement.
#' @param conf.level Desired coverage of the confidence interval for the NaRnEA Proportional Enrichment Score. Default of 0.95.
#' @param boot.num Number of times to sample the gene set members with replacement when estimating the confidence interval for the NaRnEA Proportional Enrichment Score. Default of 100. 
#' @param sig.figs Number of significant figures to display when plotting numerical values. Default of 4.
#' @return This function returns the NaRnEA Proportional Enrichment Score (PES) with a confidence interval (if desired), the NaRnEA Normalized Enrichment Score (NES), the NaRnEA two-sided p-value of the NaRnEA, and a figure made with ggplot2 that visualizes the enrichment of the gene set members in the nonparametrically transformed differential gene expression signature. 
#' @export

EnrichmentPlot <- function(signature, association.weight, association.mode, minimum.size = 30, seed = 1, title.text = "Nonparametric analytical Rank-based Enrichment Analysis", x.label.text = "Nonparametric Gene Expression Signature", line.width = 1, lower.color = "blue1", middle.color = "grey50", upper.color = "red1", conf.int = TRUE, conf.level = 0.95, boot.num = 100, sig.figs = 4){
	
	# set the seed if ggplot2 can be loaded and return a warning if not
	cur.check.value <- suppressWarnings(require("ggplot2", quietly = TRUE))
	
	if(cur.check.value){
		library(ggplot2)
		set.seed(seed)
	} else {
		stop("... install package ggplot2 from CRAN for plotting ...")
	}
	
	# set the input values to local variables
	cur.sig <- signature
	cur.rc.values <- association.weight
	cur.mor.values <- association.mode
	
	# check that the Association Weight values are positive
	if(!prod(cur.rc.values > 0)){
		stop("... association.weight values must be positive ...")
	} 
	
	# check that the Association Mode values are between (-1) and (1)
	if(!prod(abs(cur.mor.values) <= 1)){
		stop("... association.mode values must be between (-1) and (1) ...")
	}
	
	# correct and Association Mode values which are exactly equal to 0 or exactly equal to 1 
	cur.mor.values[which(abs(cur.mor.values) == 1)] <- 0.999*sign(cur.mor.values[which(abs(cur.mor.values) == 1)])
	cur.mor.values[which(cur.mor.values == 0)] <- sample(c(0.001, -0.001), size = length(which(cur.mor.values == 0)), replace = TRUE, prob = c(mean(cur.mor.values[which(cur.mor.values != 0)] > 0), mean(cur.mor.values[which(cur.mor.values != 0)] < 0)))
	
	# check that the names for the Association Weight and Association Mode match and return an error message if they do not
	if(identical(names(cur.rc.values),names(cur.mor.values))){
		cur.target.values <- names(cur.rc.values)
	} else {
		stop("... association.weight and association.mode must have the same names ...")
	}
	
	# compute the nonparametric signature transformation if no missing values or repeated entries are present; if missing values or repeated entries are present in the signature, return an error message
	if(anyNA(cur.sig)){
		stop("... signature cannot contain missing values ...")
	} else if(!prod(!duplicated(names(cur.sig)))){
		stop("... signature cannot contain duplicate entries ...")
	} else {
		set.seed(seed)
		cur.sig[which(cur.sig == 0)] <- runif(n = sum(cur.sig == 0), min = 0, max = min(abs(cur.sig[which(cur.sig != 0)])))*sample(c(-1,1), size = sum(cur.sig == 0), replace = TRUE, prob = c(mean(cur.sig[which(cur.sig != 0)] < 0),mean(cur.sig[which(cur.sig != 0)] > 0)))
		cur.np.sig <- rank(x = abs(as.numeric(cur.sig)), ties.method = "random")*sign(as.numeric(cur.sig))
		names(cur.np.sig) <- names(cur.sig)
	}
	
	# return an error message if the gene set is below the minimum size; if not, compute the indices for each gene set member
	if(length(cur.target.values) < minimum.size){
		stop("... gene set is too small ...")
	} else {
		cur.idx.values <- match(cur.target.values,names(cur.np.sig))
	}
	
	# compute the enrichment with NaRnEA
	cur.narnea.res.list <- NaRnEA(signature = signature, association.weight = association.weight, association.mode = association.mode, ledge = FALSE, minimum.size = minimum.size, seed = seed)
	names(cur.narnea.res.list) <- c("pes.value","nes.value","log.p.value")
	cur.narnea.res.list$log.p.value <- (pnorm(q = abs(cur.narnea.res.list$nes.value), lower.tail = FALSE, log.p = TRUE) + log(2))
	cur.narnea.res.list$p.value <- exp(cur.narnea.res.list$log.p.value)
	if(cur.narnea.res.list$p.value > 0){
		cur.narnea.res.list$display.p.value <- as.character(signif(cur.narnea.res.list$p.value, sig.figs))
	} else {
		cur.narnea.res.list$display.p.value <- paste("1e", ceiling(cur.narnea.res.list$log.p.value/log(10)), sep = "")
	}
	cur.narnea.res.list$display.p.value <- strsplit(x = cur.narnea.res.list$display.p.value, split = "e", fixed = TRUE)[[1]]
	if(length(cur.narnea.res.list$display.p.value) > 1){
		cur.narnea.res.list$display.p.value <- paste(substr(cur.narnea.res.list$display.p.value[1], start = 1, stop = min(nchar(cur.narnea.res.list$display.p.value[1]),sig.figs)), cur.narnea.res.list$display.p.value[2], sep = "e")
	}
	
	# if desired, estimate the confidence interval for the NaRnEA Proportional Enrichment Score using bootstrapping and a Fisher Transformation
	if(conf.int){
		cur.boot.idx.list <- lapply(1:boot.num,function(i){
			y <- sample(1:length(cur.target.values), replace = TRUE, size = length(cur.target.values))
			return(y)
		})
		cur.boot.pes.values <- sapply(cur.boot.idx.list, function(sub.idx.values){
			sub.rc.values <- cur.rc.values[sub.idx.values]
			sub.mor.values <- cur.mor.values[sub.idx.values]
			sub.pes.value <- NaRnEA(signature = signature, association.weight = sub.rc.values, association.mode = sub.mor.values, ledge = FALSE, minimum.size = minimum.size, seed = seed)$pes
			return(sub.pes.value)
		})
		cur.pes.lower.bound <- tanh((atanh(cur.narnea.res.list$pes.value) + qt(p = ((1 - conf.level)/2), lower.tail = TRUE, df = (boot.num - 1), log.p = FALSE)*sd(atanh(cur.boot.pes.values))))
		cur.pes.upper.bound <- tanh((atanh(cur.narnea.res.list$pes.value) + qt(p = ((1 - conf.level)/2), lower.tail = FALSE, df = (boot.num - 1), log.p = FALSE)*sd(atanh(cur.boot.pes.values))))
		cur.pes.conf.int <- c(cur.pes.lower.bound, cur.pes.upper.bound)
		cur.pes.conf.int[1] <- max(c(cur.pes.conf.int[1],-1))
		cur.pes.conf.int[2] <- min(c(cur.pes.conf.int[2],1))
		cur.pes.conf.int[which(is.na(cur.pes.conf.int))] <- cur.narnea.res.list$pes.value
	}
	
	# plot the targets in the nonparametrically transformed gene expression signature
	cur.plot.data <- data.frame(x.values = (cur.np.sig[cur.idx.values]/length(cur.np.sig)), mor.values = cur.mor.values, rc.values = cur.rc.values)
	
	cur.plot.data$y.start.values <- 0.5
	cur.plot.data$y.stop.values <- cur.plot.data$y.start.values + 0.5*sign(cur.plot.data$mor.values)
	
	x.breaks <- seq(from = -1, to = 1, by = 0.25)
	y.breaks <- seq(from = 0, to = 1, by = .1)
	
	cur.title <- title.text
	if(conf.int){
		cur.subtitle <- paste("PES = ", signif(cur.narnea.res.list$pes.value, sig.figs), " [",signif(cur.pes.conf.int[1], sig.figs),",",signif(cur.pes.conf.int[2], sig.figs),"]", " : ", "NES = ", signif(cur.narnea.res.list$nes.value, sig.figs), " : ", "p = ", cur.narnea.res.list$display.p.value, sep = "")
	} else {
		cur.subtitle <- paste("PES = ", signif(cur.narnea.res.list$pes.value, sig.figs), " : ", "NES = ", signif(cur.narnea.res.list$nes.value, sig.figs), " : ", "p = ", cur.narnea.res.list$display.p.value, sep = "")
	}
	cur.x.lab <- x.label.text
	
	cur.plot <- ggplot(cur.plot.data) + theme_bw() + geom_segment(x = 0, xend = 0, y = 0, yend = 1, colour = "black", lwd = 1, lineend = "square") + geom_segment(x = -1, xend = -1, y = 0, yend = 1, colour = "black", lwd = 1, lineend = "square") + geom_segment(x = 1, xend = 1, y = 0, yend = 1, colour = "black", lwd = 1, lineend = "square") + geom_segment(x = -1, xend = 1, y = 0, yend = 0, colour = "black", lwd = 1, lineend = "square") + geom_segment(x = -1, xend = 1, y = 1, yend = 1, colour = "black", lwd = 1, lineend = "square") + geom_segment(aes(x = x.values, xend = x.values, y = y.start.values, yend = y.stop.values, colour = mor.values, alpha = rc.values), lwd = line.width) + scale_colour_gradient2(low = lower.color, mid = middle.color, high = upper.color, midpoint = 0) + scale_alpha_continuous(limits = c((1/nrow(cur.plot.data)),1)) + guides(colour = "none", alpha = "none") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, colour = "black", size = 25), plot.subtitle = element_text(hjust = 0.5, colour = "black", size = 22), axis.title.x = element_text(hjust = 0.5, colour = "black", size = 20), axis.text.x = element_text(hjust = 0.5, colour = "black", size = 18), plot.margin = ggplot2::margin(t = 10, r = 25, b = 5, l = 25, unit = "pt")) + scale_x_continuous(limits = range(x.breaks), breaks = x.breaks, expand = c(0,0)) + scale_y_continuous(limits = range(y.breaks), breaks = y.breaks, expand = c(0,0)) + labs(title = cur.title, subtitle = cur.subtitle, x = cur.x.lab)
	
	return(cur.plot)
	
}






