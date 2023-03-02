# NaRnEA

## NaRnEA (Nonparametric analytical Rank-based Enrichment Analysis)

Nonparametric analytical Rank-based Enrichment Analysis (`NaRnEA`) is a novel statistical method created by Griffin et al. for the purpose of performing powerful and accurate gene set analysis. To read more about the theory behind `NaRnEA` and how it compares with Gene Set Enrichment Analysis (`GSEA`) and analytical Rank-based Enrichment Analysis (`aREA`), check out the preprint of the `NaRnEA` manuscript on bioRxiv ([link](https://www.biorxiv.org/content/10.1101/2021.10.02.462873v2.full)).

## Installing the NaRnEA R package

To install the `NaRnEA` R package, download the zipped GitHub codebase, unzip it, and install locally with the `devtools` package using the code below (if the unzipped directory is named `NaRnEA-main`, change this to `NaRnEA` to ensure subsequent code functionality). The `NaRnEA` GitHub codebase also includes proteomic data and code for running `ARACNe3`, the newest implementation of the Algorithm for the Reconstruction of Accurate Cellular Networks, so that users can replicate the analysis from the `NaRnEA` manuscript. The code is freely available for academic research use; for more information, see the license associated with the `NaRnEA` GitHub codebase. To run `ARACNe3` users will need to install the latest version of JDK.

```{r}
devtools::install(pkg = "NaRnEA", build_vignettes = TRUE)
```

## Using NaRnEA and ARACNe3
`NaRnEA` performs gene set analysis using two objects; a gene expression signature, and a list of gene sets parametrized with directionality (association mode) and confidence / importance (association weight) values. Below is a simple example of how to run `NaRnEA` if you have these two objects already generated:

```{r}
library(NaRnEA)
library(ggplot2)
data("LUAD_ges")
data("LUAD_network")

## run narnea
luad.narnea.res <- lapply(LUAD_network, function(gene.set, cur.ges) {
  sub.aw.values <- as.numeric(gene.set$aw.values)
  names(sub.aw.values) <- as.character(gene.set$target.values)
  
  sub.am.values <- as.numeric(gene.set$am.values)
  names(sub.am.values) <- as.character(gene.set$target.values)
  
  sub.narnea.res <- NaRnEA(signature = cur.ges, 
                           association.weight = sub.aw.values, 
                           association.mode = sub.am.values, 
                           ledge = FALSE)
  
  return(sub.narnea.res)
}, cur.ges = LUAD_ges)

## visualize NaRnEA results with a volcano plot
# make plot data
plot.data <- as.data.frame(do.call("rbind", lapply(luad.narnea.res, function(x){
  y <- c("pes.values" = x$pes, "nes.values" = x$nes)
  return(y)
})))
plot.data$p.values <- 2*pnorm(q = abs(plot.data$nes.values), lower.tail = FALSE)
plot.data$fdr.values <- p.adjust(p = plot.data$p.values, method = "BH")
plot.data$dir.values <- "NS"
plot.data[which(plot.data$pes.values > 0 & plot.data$fdr.values < 0.05),"dir.values"] <- "UP"
plot.data[which(plot.data$pes.values < 0 & plot.data$fdr.values < 0.05),"dir.values"] <- "DOWN"
plot.data$dir.values <- factor(plot.data$dir.values, levels = c("DOWN", "NS", "UP"))
# set plot parameters
colour.value.vec <- c("DOWN" = "blue1", "NS" = "grey50", "UP" = "red1")
x.breaks <- seq(from = -1, to = 1, by = 0.2)
y.breaks <- seq(from = 0, to = 45, by = 5)
cur.title <- "NaRnEA Results - LUAD"
cur.subtitle <- paste("DOWN = ", signif((100*mean(plot.data$dir.values == "DOWN")), 4), "% (",
                      sum(plot.data$dir.values == "DOWN")," / ",nrow(plot.data),")", 
                      " : ", "UP = ", signif((100*mean(plot.data$dir.values == "UP")), 4), "% (",
                      sum(plot.data$dir.values == "UP")," / ",nrow(plot.data),")", sep = "")
cur.x.lab <- "Proportional Enrichment Score"
cur.y.lab <- "| Normalized Enrichment Score |"
cur.color.lab <- "Enrichment (FDR < 0.05)"
# generate plot
ggplot(plot.data) + theme_bw() + 
  geom_hline(colour = "black", lwd = 1, yintercept = 0) + 
  geom_vline(colour = "black", lwd = 1, xintercept = 0) + 
  geom_point(aes(x = pes.values, y = abs(nes.values), colour = dir.values), 
             size = 3, shape = 21, fill = "white", alpha = 1, stroke = 1) + 
  scale_colour_manual(values = colour.value.vec) + 
  scale_x_continuous(limits = range(x.breaks), breaks = x.breaks) + 
  scale_y_continuous(limits = range(y.breaks), breaks = y.breaks) + 
  labs(title = cur.title, subtitle = cur.subtitle, x = cur.x.lab, y = cur.y.lab, colour = cur.color.lab) + 
  theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 30), 
        plot.subtitle = element_text(hjust = 0.5, colour = "black", size = 23), 
        axis.title = element_text(hjust = 0.5, colour = "black", size = 25), 
        axis.text.x = element_text(hjust = 0.5, colour = "black", size = 22), 
        axis.text.y = element_text(hjust = 1, colour = "black", size = 22), 
        legend.title = element_text(hjust = 0, colour = "black", size = 20), 
        legend.text = element_text(hjust = 0, colour = "black", size = 18), 
        legend.position = "bottom")	

## visualize a positively enriched target (FOXM1) from the NaRnEA results using EnrichmentPlot
pos.reg <- 'g_2305_'
gene.set <- LUAD_network[[pos.reg]]
# filter am and aw values
sub.aw.values <- as.numeric(gene.set$aw.values)
names(sub.aw.values) <- as.character(gene.set$target.values)
sub.am.values <- as.numeric(gene.set$am.values)
names(sub.am.values) <- as.character(gene.set$target.values)
# generate plot
EnrichmentPlot(LUAD_ges, sub.aw.values, sub.am.values, title.text = 'FOXM1 NaRnEA Enrichment')

## visualize a negatively enriched target (TCF21) from the NaRnEA results using EnrichmentPlot
pos.reg <- 'g_6943_'
gene.set <- LUAD_network[[pos.reg]]
# filter am and aw values
sub.aw.values <- as.numeric(gene.set$aw.values)
names(sub.aw.values) <- as.character(gene.set$target.values)
sub.am.values <- as.numeric(gene.set$am.values)
names(sub.am.values) <- as.character(gene.set$target.values)
# generate plot
EnrichmentPlot(LUAD_ges, sub.aw.values, sub.am.values, title.text = 'TCF21 NaRnEA Enrichment')
```

For a more full walkthrough that covers both generating the gene expression signature and gene set objects (using ARACNe3), please see the full vignettes included in the `NaRnEA` package. These vignettes cover all analysis performed in the `NaRnEA` manuscript, including the reverse-engineering of transcriptional regulatory networks from gene expression data in The Cancer Genome Atlas (TCGA) for lung adenocarcinoma (LUAD), colon adenocarcinoma (COAD), and head-neck squamous cell carcinoma (HSNC) cancer subtypes with `ARACNe3`. 
To learn more about how to use `NaRnEA` for performing powerful and accurate gene set analysis, check out the `NaRnEA` vignettes with the following code.

```{r}
library(NaRnEA)
browseVignettes(package = "NaRnEA")
```
