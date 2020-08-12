# NaRnEA
## NaRnEA (Non-parametric analytical rank-based enrichment analysis)

Non-parametric analytical Rank-based Enrichment Analysis (NaRnEA) enables users to perform highly accurate gene set analysis in a non-parametric, frequentist manner.

## Simulated Gene Set Enrichment Walkthrough 

install the NaRnEA package from github
```{r}
library(devtools)
install_github(repo = "califano-lab/NaRnEA", auth_token = "25fa79a6847c853f310a01dcfe173091ccb3400a", force = TRUE)
```

load the NaRnEA package and ggplot2 for plotting
```{r}
library(NaRnEA)
library(ggplot2)
```

set the simulation parameters
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

set the seed for the simulation 
```{r}
set.seed(sim.seed)
```

simulate the gene expression signature values for all genes in the different samples under the null distribution
```{r}
cur.ges.mat <- matrix(data = NA, nrow = gene.num, ncol = sample.num)
cur.ges.mat <- apply(cur.ges.mat,2,function(x){
	y <- runif(n = length(x), min = null.ges.min, max = null.ges.max)
	return(y)
})
colnames(cur.ges.mat) <- paste("s",1:ncol(cur.ges.mat),"",sep = "_")
rownames(cur.ges.mat) <- paste("g",1:nrow(cur.ges.mat),"",sep = "_")
```

randomly select some genes and construct the gene sets which do not overlap
```{r}
gene.set.1.targets <- sample(rownames(cur.ges.mat), size = gene.set.size)
gene.set.2.targets <- sample(setdiff(rownames(cur.ges.mat),gene.set.1.targets), size = gene.set.size)

gene.set.1 <- list(tfmode = runif(n = gene.set.size, min = gene.set.MoR.min, max = gene.set.MoR.max), likelihood = runif(n = gene.set.size, min = gene.set.RC.min, max = gene.set.RC.max))
names(gene.set.1$tfmode) <- gene.set.1.targets

gene.set.2 <- list(tfmode = runif(n = gene.set.size, min = gene.set.MoR.min, max = gene.set.MoR.max), likelihood = runif(n = gene.set.size, min = gene.set.RC.min, max = gene.set.RC.max))
names(gene.set.2$tfmode) <- gene.set.2.targets
```

verify the gene sets don't overlap
```{r}
intersect(names(gene.set.1$tfmode),names(gene.set.2$tfmode))
```

modulate the gene expression signature values for the members of the first gene set
```{r}
cur.ges.mat <- apply(cur.ges.mat,2,function(x){
	y <- x
	y[match(names(gene.set.1$tfmode),names(x))] <- runif(n = length(gene.set.1$tfmode), min = alt.ges.min, max = alt.ges.max)*sign(as.numeric(gene.set.1$tfmode))
	return(y)
})
```

compute the enrichment of the gene sets in the gene expression signatures using NaRnEA
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


 
