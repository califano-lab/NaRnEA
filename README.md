# NaRnEA
## NaRnEA (Non-parameetric analytical rank-based enrichment analysis)

Non-parametric analytical Rank-based Enrichment Analysis (NaRnEA) enables users to perform highly accurate gene set analysis in a non-parametric, frequentist manner.

## Walkthrough 

Use the following code to download the NaRnEA R package

```{r}
devtools::install_github("califano-lab/NaRnEA")
```

Load the NaRnEA package
```{r}
library(NaRnEA)
```

Now we'll set the global parameters for the simulation
```{r}
sim.seed <- 11
gene.num <- 20000
null.ges.min <- -10
null.ges.max <- -10
alt.ges.min <- -10
alt.ges.max <- 10
gene.set.size <- 100
gene.set.RC.min <- 0
gene.set.RC.max <- 1
gene.set.MoR.min <- -1
gene.set.MoR.max <- 1
```


 
