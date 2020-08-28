# NaRnEA

## NaRnEA (Non-parametric analytical rank-based enrichment analysis)

Non-parametric analytical Rank-based Enrichment Analysis (NaRnEA) enables users to perform accurate and powerful gene set analysis in a non-parametric, frequentist manner.

## Installing NaRnEA

NaRnEA is distributed free for scientific use as an R package; use the following code to download and install NaRnEA.
```{r}
devtools::install_github(repo = "califano-lab/NaRnEA")
```

## Simulated Gene Set Enrichment Analysis

Load the NaRnEA package and the ggplot2 package (for improved plotting functionality)
```{r}
library(NaRnEA)
library(ggplot2)
```

Set the parameters for the simulation. Here we will be simulating 5,000 gene expression signatures, each with 20,000 genes. The gene expression signature values for each gene will be independently and identically distributed from a uniform distribution on the interval [-10,10]. NaRnEA is a non-parametric method and therefore does not assume anything about the gene expression signature, but using a simple gene expression signature for these simulations increases our ability to intuitively interpret the results. 
```{r}
sim.seed <- 1234
gene.num <- 20000
sample.num <- 5000
sig.figs <- 3
boot.num <- 1E3
null.ges.min <- -10
null.ges.max <- 10
alt.ges.min <- 5
alt.ges.max <- 10
gene.set.size <- 100
gene.set.RC.min <- 0
gene.set.RC.max <- 1
gene.set.MoR.min <- -1
gene.set.MoR.max <- 1
```

Set the seed for the simulation; this will ensure that the results of the simulation are reproducible.
```{r}
set.seed(sim.seed)
```

Now we'll simulate the gene expression signature values and store them in a matrix with genes in the rows and samples in the columns.
```{r}
cur.ges.mat <- matrix(data = NA, nrow = gene.num, ncol = sample.num)
cur.ges.mat <- apply(cur.ges.mat,2,function(x){
	y <- runif(n = length(x), min = null.ges.min, max = null.ges.max)
	return(y)
})
colnames(cur.ges.mat) <- paste("s",1:ncol(cur.ges.mat),"",sep = "_")
rownames(cur.ges.mat) <- paste("g",1:nrow(cur.ges.mat),"",sep = "_")
```

Next we will create two gene sets, each with the prespecified number of members. For this simulation the gene sets will be constructed so that they do not overlap. The Regulation Confidence (`r likelihood`) values will be randomly drawn from a uniform distribution bounded by the prespecificed `r gene.set.RC.min` and `r gene.set.RC.max` parameters. The Mode of Regulation (`r tfmode`) values will be randomly drawn from a uniform distribution bounded by the prespecificed `r gene.set.MoR.min` and `r gene.set.MoR.max` parameters. 
```{r}
gene.set.1.targets <- sample(rownames(cur.ges.mat), size = gene.set.size)
gene.set.2.targets <- sample(setdiff(rownames(cur.ges.mat),gene.set.1.targets), size = gene.set.size)

gene.set.1 <- list(tfmode = runif(n = gene.set.size, min = gene.set.MoR.min, max = gene.set.MoR.max), likelihood = runif(n = gene.set.size, min = gene.set.RC.min, max = gene.set.RC.max))
names(gene.set.1$tfmode) <- gene.set.1.targets

gene.set.2 <- list(tfmode = runif(n = gene.set.size, min = gene.set.MoR.min, max = gene.set.MoR.max), likelihood = runif(n = gene.set.size, min = gene.set.RC.min, max = gene.set.RC.max))
names(gene.set.2$tfmode) <- gene.set.2.targets
```

We can verify quickly whether or not the gene sets overlap.
```{r}
intersect(names(gene.set.1$tfmode),names(gene.set.2$tfmode))
```

As of now neither gene set is enriched in any gene expression signature; all gene set members are independently and identically distributed in the gene expression signatures. We will now change that by resimulating the gene expression signature values for the members of the first gene set based on the `r alt.ges.min` and `r alt.ges.max` parameters and the sign of each target's Mode of Regulation.
```{r}
cur.ges.mat <- apply(cur.ges.mat,2,function(x){
	y <- x
	y[match(names(gene.set.1$tfmode),names(x))] <- runif(n = length(gene.set.1$tfmode), min = alt.ges.min, max = alt.ges.max)*sign(as.numeric(gene.set.1$tfmode))
	return(y)
})
```

We can now compute the enrichment of the gene sets in each of the gene expression signatures with NaRnEA. For now, we'll do this without leading edge analysis; setting `r leading.edge = TRUE` allows for post-hoc leading edge analysis but increases computational time. 
```{r}
gene.set.1.enrichment.results <- apply(cur.ges.mat,2,function(cur.ges,cur.gene.set){
	y <- unlist(NaRnEA(ges = cur.ges, regulon = cur.gene.set, seed = 1, leading.edge = FALSE))
	return(y)
}, cur.gene.set = gene.set.1)

gene.set.2.enrichment.results <- apply(cur.ges.mat,2,function(cur.ges,cur.gene.set){
	y <- unlist(NaRnEA(ges = cur.ges, regulon = cur.gene.set, seed = 1, leading.edge = FALSE))
	return(y)
}, cur.gene.set = gene.set.2)
```

reformat the results as data frames and compute two-sided p-values, false discovery rate values, and family wise error rate values
```{r}
gene.set.1.enrichment.results <- as.data.frame(t(gene.set.1.enrichment.results))
gene.set.1.enrichment.results$p.values <- 2*pnorm(q = abs(gene.set.1.enrichment.results$nes), lower.tail = FALSE)
gene.set.1.enrichment.results$fdr.values <- p.adjust(p = gene.set.1.enrichment.results$p.values, method = "BH")
gene.set.1.enrichment.results$fwer.values <- p.adjust(p = gene.set.1.enrichment.results$p.values, method = "bonferroni")

gene.set.2.enrichment.results <- as.data.frame(t(gene.set.2.enrichment.results))
gene.set.2.enrichment.results$p.values <- 2*pnorm(q = abs(gene.set.2.enrichment.results$nes), lower.tail = FALSE)
gene.set.2.enrichment.results$fdr.values <- p.adjust(p = gene.set.2.enrichment.results$p.values, method = "BH")
gene.set.2.enrichment.results$fwer.values <- p.adjust(p = gene.set.2.enrichment.results$p.values, method = "bonferroni")
```

visualize the results 
```{r}
plot.data <- data.frame(pes.values = c(gene.set.1.enrichment.results$pes, gene.set.2.enrichment.results$pes), nes.values = c(gene.set.1.enrichment.results$nes, gene.set.2.enrichment.results$nes), p.values = c(gene.set.1.enrichment.results$p.values, gene.set.2.enrichment.results$p.values), group.values = c(rep("GS1", times = nrow(gene.set.1.enrichment.results)), rep("GS2", times = nrow(gene.set.1.enrichment.results))))

plot.data$group.values <- factor(plot.data$group.values, levels = c("GS1","GS2"), labels = c("Gene Set 1", "Gene Set 2"))

colour.value.vec <- c("Gene Set 1" = "red1", "Gene Set 2" = "blue1")

x.step <- .1
x.min <- -1
x.max <- 1
x.breaks <- seq(from = x.min, to = x.max, by = x.step)

bin.adjust <- .1
bin.breaks <- seq(from = x.min, to = x.max, by = (bin.adjust*x.step))

gene.set.1.binom.res <- binom.test(x = sum(gene.set.1.enrichment.results$p.values < .05), n = nrow(gene.set.1.enrichment.results), p = .05, alternative = "greater")
gene.set.2.binom.res <- binom.test(x = sum(gene.set.2.enrichment.results$p.values < .05), n = nrow(gene.set.2.enrichment.results), p = .05, alternative = "greater")

cur.title <- "Simulated Gene Set Enrichment"
cur.subtitle <- paste("Genes = ", gene.num, " : ", "Samples = ", sample.num, " : ", "Gene Set Size = ", gene.set.size, "\n", "Null GES ~ U[",null.ges.min,",",null.ges.max,"]", " : ", "Alt GES ~ U[",alt.ges.min,",",alt.ges.max,"]", "\n", "Gene Expression Signatures with Statistically Significant Enrichment", "\n", "Gene Set 1 = ",signif(100*as.numeric(gene.set.1.binom.res$estimate),sig.figs), "% (p = ",signif(as.numeric(gene.set.1.binom.res$p.value),sig.figs),")", " : ", "Gene Set 2 = ",signif(100*as.numeric(gene.set.2.binom.res$estimate),sig.figs), "% (p = ",signif(as.numeric(gene.set.2.binom.res$p.value),sig.figs),")",sep = "")
cur.x.lab <- "NaRnEA Proportional Enrichment Scores (PES)"
cur.y.lab <- "Probability Density"

ggplot(plot.data) + facet_wrap(~group.values, ncol = 1, scales = "free") + geom_vline(colour = "black", lwd = .5, xintercept = c(0,range(x.breaks))) + geom_histogram(aes(x = pes.values, y = ..density.., colour = group.values), lwd = .5, fill = "white", breaks = bin.breaks, closed = "left") + scale_colour_manual(values = colour.value.vec) + scale_x_continuous(limits = range(x.breaks), breaks = x.breaks) + geom_hline(colour = "black", lwd = .5, yintercept = 0) + labs(title = cur.title, subtitle = cur.subtitle, x = cur.x.lab, y = cur.y.lab) + theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 25), plot.subtitle = element_text(hjust = 0.5, colour = "black", size = 22), axis.title = element_text(hjust = 0.5, colour = "black", size = 20), strip.text = element_text(hjust = 0.5, colour = "black", size = 20), axis.text.x = element_text(hjust = 0.5, colour = "black", size = 18), axis.text.y = element_text(hjust = 1, colour = "black", size = 18)) + guides(colour = FALSE)
```


 
