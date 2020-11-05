#' Enrichment Analysis Plotting function
#' 
#' @param cur.ges Gene expression signature (named numeric vector)
#' @param cur.regul Gene set parameterized with Regulation Confidence (likelihood) and Mode of Regulation (tfmode)
#' @param cur.minsize Minimum number of gene set members that are present in the gene expression signature. Default of 30.
#' @param cur.seed Random number generator seed to ensure ties in the gene expression signature are broken in a reproducible manner. Default of 1.
#' @param cur.title The text for the main title.
#' @param cur.subtitle The text for the subtitle. Default is NULL. If non-NULL, this will replace the output NaRnEA PES, NES, and two-sided p-value for the enrichment of the gene set in the gene expression signature.
#' @param cur.x.lab The text for the x-axis label.
#' @param cur.y.lab The text for the y-axis label.
#' @param title.font.size The font size for the main title.
#' @param subtitle.font.size The font size for the subtitle.
#' @param axis.title.font.size The font size for the x-axis label and the y-axis label.
#' @param axis.text.font.size The font size for the text of the x-axis and y-axis.
#' @param positive.color The color for the plot elements corresponding to gene set members with a positive Mode of Regulation (tfmode).
#' @param middle.color The middle color for the color gradient between the positive color and the negative color.
#' @param negative.color The color for the plot elements corresponding to gene set members with a negative Mode of Regulation (tfmode)
#' @param enrich.lwd The linewidth for the running sum statistic visualization.
#' @param mid.bar.lwd The linewidth for the horizontal bar denoting the expected value of the running sum statistic.
#' @param black.bar.lwd The linewidth for the vertical bar denoting where the gene expression signature switches from negative values to positive values. 
#' @param sig.lwd The linewidth for the vertical lines denoting the location of each gene set member in the gene expression signature.
#' @return This function returns the NaRnEA PES, NES, and two-sided p-value for the enrichment of the gene set in the gene expression signature as well as a figure made with ggplot2 that visualizes the enrichment of the gene set in the gene expression signature using the two-sample Kolmogorov Smirnov running sum statistic (if cur.score = 0) or the GSEA running sum statistic (if cur.score = 1).
#' @export
EnrichmentAnalysisPlot <- function(cur.ges, cur.regul, cur.minsize = 30, cur.score = 0, cur.seed = 1, non.parametric = TRUE, cur.title = "Non-parametric Analytical Rank-based Enrichment Analysis", cur.subtitle = NULL, cur.x.lab = "Gene Expression Signature", cur.y.lab = "KS Running Sum Statistic", title.font.size = 25, subtitle.font.size = 22, axis.title.font.size = 20, axis.text.font.size = 18, positive.color = "red1", middle.color = "grey40", negative.color = "blue1", enrich.lwd = 1, mid.bar.lwd = 1, black.bar.lwd = .75, sig.lwd = .85){
	
	# set the seed
	set.seed(cur.seed)
	
	# correct the gene expression signature
	cur.ges <- cur.ges[is.finite(cur.ges)]
	cur.ges <- cur.ges[which(cur.ges != 0)]
	
	# correct the regulon as necessary
	if(prod(abs(cur.regul$tfmode) == 1)){
		cur.regul$tfmode <- cur.regul$tfmode*(.99)
	}
	
	if(anyNA(match(names(cur.regul$tfmode),names(cur.ges)))){
		keep.index <- which(names(cur.regul$tfmode)%in%names(cur.ges))
		cur.regul$tfmode <- cur.regul$tfmode[keep.index]
		cur.regul$likelihood <- cur.regul$likelihood[keep.index]
	}
	
	
	# set the geneset for the KS enrichment score calculations
	cur.geneset <- cur.regul$tfmode
	
	# compute the NaRnEA normalized enrichment score, proportional enrichment score, and statistical significance
	cur.narnea.res <- NaRnEA(ges = cur.ges, regulon = cur.regul, minsize = cur.minsize, leading.edge = FALSE)
	cur.nes <- cur.narnea.res$nes
	cur.pes <- cur.narnea.res$pes
	
	narnea.nes.value <- signif(cur.nes,3)
	narnea.pes.value <- signif(min(cur.pes,1),3)
	narnea.log10.pvalue <- (log(2) + pnorm(q = abs(cur.nes), lower.tail = FALSE, log.p = TRUE))/log(10)
	narnea.pvalue.power <- floor(narnea.log10.pvalue)
	narnea.pvalue.num <- signif((10^(narnea.log10.pvalue - floor(narnea.log10.pvalue))),3)
	if(is.null(cur.subtitle)){
		if(2*pnorm(q = abs(cur.nes), lower.tail = FALSE, log.p = FALSE) >= 0.01){
			cur.subtitle <- paste("NaRnEA"," : ","PES = ",narnea.pes.value," : ","NES = ",narnea.nes.value," : ", "p = ",signif(2*pnorm(q = abs(cur.nes), lower.tail = FALSE, log.p = FALSE),3),sep = "")
		} else {
			cur.subtitle <- paste("NaRnEA"," : ","PES = ",narnea.pes.value," : ","NES = ",narnea.nes.value," : ","p = ",narnea.pvalue.num,"e",narnea.pvalue.power,sep = "")
		}
	}
	
	# perform a non-parametric transformation of the gene expression signature for the KS enrichment analysis if desired
	cur.sig <- cur.ges
	if(non.parametric){
			cur.sig <- sign(cur.sig)*rank(abs(cur.sig),ties.method = "random")
	}
	
	# compute the KS enrichment score for the positive targets
	if(max(cur.geneset) > 0){
		cur.up.sig <- sort(cur.sig,decreasing = TRUE)
		cur.up.geneset <- cur.geneset[which(cur.geneset > 0)]
		cur.up.index <- match(names(cur.up.geneset),names(cur.up.sig))
		cur.up.es.values <- -1*sign(abs(cur.up.sig))/(length(cur.up.sig) - length(cur.up.geneset))
		cur.up.es.values[cur.up.index] <- (abs(cur.up.sig[cur.up.index])^(cur.score))/sum((abs(cur.up.sig[cur.up.index])^(cur.score)))
		cur.up.es.values <- cumsum(cur.up.es.values)
		cur.up.es <- as.numeric(cur.up.es.values[which(abs(cur.up.es.values) == max(abs(cur.up.es.values)))])
		cur.up.es <- cur.up.es[1]
		cur.up.es.values[1] <- 0
		cur.up.es.values[length(cur.up.es.values)] <- 0
	} else {
		cur.up.es <- 0
		cur.up.es.values <- rep(0,times = length(cur.sig))
	}
	
	# compute the KS enrichment score for the negative targets
	if(min(cur.geneset) < 0){
		cur.down.sig <- sort(cur.sig,decreasing = FALSE)
		cur.down.geneset <- cur.geneset[which(cur.geneset < 0)]
		cur.down.index <- match(names(cur.down.geneset),names(cur.down.sig))
		cur.down.es.values <- -1*sign(abs(cur.down.sig))/(length(cur.down.sig) - length(cur.down.geneset))
		cur.down.es.values[cur.down.index] <- (abs(cur.down.sig[cur.down.index])^(cur.score))/sum(abs(cur.down.sig[cur.down.index])^(cur.score))
		cur.down.es.values <- -1*cumsum(cur.down.es.values)
		cur.down.es <- as.numeric(cur.down.es.values[which(abs(cur.down.es.values) == max(abs(cur.down.es.values)))])
		cur.down.es <- cur.down.es[1]
		cur.down.es.values[1] <- 0
		cur.down.es.values[length(cur.down.es.values)] <- 0
	} else {
		cur.down.es <- 0
		cur.down.es.values <- rep(0,times = length(cur.sig))
	}
	
	
	# plot the enrichment
	
	enrich.plot.data <- data.frame(x.values = rep(1:length(cur.sig),times = 2), y.values = c(rev(cur.up.es.values),(cur.down.es.values)),type = c(rep("UP",times = length(cur.sig)),rep("DOWN",times = length(cur.sig))), gene = c(names(rev(cur.up.es.values)),names(cur.down.es.values)))
	enrich.plot.data$type = factor(enrich.plot.data$type, levels = c("UP","DOWN"))
	enrich.plot.data$sig.values <- cur.sig[match(enrich.plot.data$gene,names(cur.sig))]
	enrich.plot.data$x.values <- enrich.plot.data$x.values/length(cur.sig)
	enrich.plot.data$MoR <- NULL
	enrich.plot.data$IC <- NULL
	enrich.plot.data[match(names(cur.regul$tfmode),enrich.plot.data$gene),"MoR"] <- as.numeric(cur.regul$tfmode)
	enrich.plot.data[match(names(cur.regul$tfmode),enrich.plot.data$gene),"IC"] <- as.numeric(cur.regul$likelihood)
	
	enrich.plot.output <- ggplot(enrich.plot.data) + theme_bw() + geom_hline(yintercept = seq(from = -1, to = 1, by = .2), colour = "grey90", lwd = .5) + geom_vline(xintercept = seq(from = 0, to = 1, by = .1), colour = "grey90", lwd = .5)
	
	if(max(cur.regul$tfmode) > 0){
		enrich.plot.output <- enrich.plot.output + geom_segment(x = enrich.plot.data[which(enrich.plot.data$y.values == cur.up.es & enrich.plot.data$type == "UP")[1],"x.values"], xend = enrich.plot.data[which(enrich.plot.data$y.values == cur.up.es & enrich.plot.data$type == "UP")[1],"x.values"], y = 0, yend = cur.up.es, colour = positive.color, lwd = enrich.lwd)
	
	}
	
	if(min(cur.regul$tfmode) < 0){
		enrich.plot.output <- enrich.plot.output + geom_segment(x = enrich.plot.data[which(enrich.plot.data$y.values == cur.down.es & enrich.plot.data$type == "DOWN")[1],"x.values"], xend = enrich.plot.data[which(enrich.plot.data$y.values == cur.down.es & enrich.plot.data$type == "DOWN")[1],"x.values"], y = 0, yend = cur.down.es, colour = negative.color, lwd = enrich.lwd)
	}
	
	if(min(cur.regul$tfmode) < 0){
		enrich.plot.output <- enrich.plot.output + geom_line(data = enrich.plot.data[which(enrich.plot.data$type == "DOWN"),],aes(x = x.values, y = y.values), colour = negative.color, lwd = enrich.lwd)
	}
	
	if(max(cur.regul$tfmode) > 0){
		enrich.plot.output <- enrich.plot.output + geom_line(data = enrich.plot.data[which(enrich.plot.data$type == "UP"),],aes(x = x.values, y = y.values), colour = positive.color, lwd = enrich.lwd)
	}
	
	enrich.plot.output <- enrich.plot.output + geom_segment(x = 0, xend = 1, y = 0, yend = 0, colour = "black", lwd = mid.bar.lwd)
	
	enrich.plot.output <- enrich.plot.output + guides(colour = FALSE, alpha = FALSE)
	
	enrich.plot.output <- enrich.plot.output + scale_y_continuous(limits = c(-1,1.26), breaks = seq(from = -1, to = 1, by = .2))
	
	enrich.plot.output <- enrich.plot.output + scale_x_continuous(limits = c(0,1), breaks = seq(from = 0, to = 1, by = .1))
	
	enrich.plot.output <- enrich.plot.output + geom_segment(colour = "black", x = 0, xend = 0, y = 1.05, yend = 1.25, lwd = black.bar.lwd)
	
	enrich.plot.output <- enrich.plot.output + geom_segment(colour = "black", x = 1, xend = 1, y = 1.05, yend = 1.25, lwd = black.bar.lwd)
	
	enrich.plot.output <- enrich.plot.output + geom_segment(colour = "black", x = mean(cur.sig < 0), xend = mean(cur.sig < 0), y = 1.05, yend = 1.25, lwd = black.bar.lwd)
	
	enrich.plot.output <- enrich.plot.output + geom_segment(x = enrich.plot.data[which(abs(enrich.plot.data$sig.values) == 1),"x.values"][1], xend = enrich.plot.data[which(abs(enrich.plot.data$sig.values) == 1),"x.values"][1], y = 1.05, yend = 1.25, colour = "black", lwd = black.bar.lwd)
	
	enrich.plot.output <- enrich.plot.output + geom_segment(data = enrich.plot.data[which(is.finite(enrich.plot.data$MoR)),],aes(x = x.values, xend = x.values, yend = 1.15, y = (1.15 + sign(MoR)/10), colour = MoR, alpha = IC), lwd = sig.lwd, inherit.aes = FALSE)
	
	enrich.plot.output <- enrich.plot.output + geom_segment(colour = "black", x = 0, xend = 1, y = 1.15, yend = 1.15, lwd = black.bar.lwd)
	
	enrich.plot.output <- enrich.plot.output + geom_segment(colour = "black", x = 0, xend = 1, y = 1.05, yend = 1.05, lwd = black.bar.lwd)
	
	enrich.plot.output <- enrich.plot.output + geom_segment(colour = "black", x = 0, xend = 1, y = 1.25, yend = 1.25, lwd = black.bar.lwd)
	
	enrich.plot.output <- enrich.plot.output + scale_colour_gradient2(low = negative.color, mid = middle.color, high = positive.color, breaks = seq(from = -1, to = 1, by = .1))
	
	enrich.plot.output <- enrich.plot.output + scale_alpha(range = c(0.5,1)) 
	
	enrich.plot.output <- enrich.plot.output + theme(plot.title = element_text(colour = "black", size = title.font.size, hjust = 0.5), plot.subtitle = element_text(colour = "black", size = subtitle.font.size, hjust = 0.5), axis.title = element_text(hjust = 0.5, colour = "black", size = axis.title.font.size), axis.text = element_text(hjust = 0.5, colour = "black", size = axis.text.font.size), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	enrich.plot.output <- enrich.plot.output + labs(title = cur.title, subtitle = cur.subtitle, x = cur.x.lab, y = cur.y.lab)
	
	print(enrich.plot.output)

}

