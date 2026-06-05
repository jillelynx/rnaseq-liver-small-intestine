---
title: "RNA-seq Differential Expression Analysis: Liver vs Small Intestine"
author: "Jill Lynch"
date: "`r Sys.Date()`"
output: html_notebook
---

## Load packages

```{r}
library(edgeR)
library(statmod)
library(limma)
library(pheatmap)
library(knitr)
```

## Load count data

```{r}
counts <- read.table("final_counts.tsv",
                     sep = "\t",
                     header = TRUE,
                     row.names = 1)
```

## Define sample groups

```{r}
group <- factor(c("Liver", "Liver", "Liver", "SI", "SI", "SI", "SI"))
```

## Create DGEList object

```{r}
y <- DGEList(counts = counts,
             genes = data.frame(Gene = rownames(counts)))

rownames(y$counts) <- rownames(y$genes) <- y$genes$Gene
```

## Filter low-expression genes and normalize

```{r}
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)
y$samples$group <- group
```

## MDS plot

```{r}
plotMDS(y,
        labels = colnames(y),
        col = c("red", "red", "red", "blue", "blue", "blue", "blue"))

title(main = "MDS Plot: Liver vs Small Intestine",
      cex.main = 1.5,
      font.main = 2)
```

## Heatmap of top 50 variable genes

```{r}
log_counts <- cpm(y, log = TRUE, prior.count = 2)

gene_variance <- apply(log_counts, 1, var)

top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:50]

heatmap_data <- log_counts[top_genes, ]

heatmap_scaled <- t(scale(t(heatmap_data)))

pheatmap(heatmap_scaled,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         main = "Top 50 Variable Genes Heatmap")
```

## Exact test analysis

```{r}
y <- estimateDisp(y)

et <- exactTest(y, pair = c("Liver", "SI"))

de <- decideTests(et)
summary(de)
```

## Biological coefficient of variation plot

```{r}
plotBCV(y)

title(main = "Biological Coefficient of Variation (BCV) Plot",
      cex.main = 1.5,
      font.main = 2)
```

## Smear plot from exact test

```{r}
detags <- rownames(y)[as.logical(de)]

plotSmear(et, de.tags = detags)
abline(h = c(-1, 1), col = "blue")

title(main = "Smear Plot: Exact Test",
      cex.main = 1.5,
      font.main = 2)
```

## Generalized linear model design

```{r}
design <- model.matrix(~ group)
rownames(design) <- colnames(y)

design
```

## GLM differential expression analysis

```{r}
y2 <- y

y2 <- estimateDisp(y2, design, robust = TRUE)

fit <- glmFit(y2, design)

lrt <- glmLRT(fit, coef = 2)

de2 <- decideTests(lrt)
summary(de2)
```

## Smear plot from GLM

```{r}
detags <- rownames(y2)[as.logical(de2)]

plotSmear(lrt, de.tags = detags)
abline(h = c(-1, 1), col = "blue")

title(main = "Smear Plot: Generalized Linear Model",
      cex.main = 1.5,
      font.main = 2)
```

## Extract full results table

```{r}
results <- topTags(lrt, n = Inf)$table

results_sorted <- results[order(results$FDR), ]

top10 <- head(results_sorted, 10)

top10$logFC <- round(top10$logFC, 2)
top10$logCPM <- round(top10$logCPM, 2)
top10$FDR <- signif(top10$FDR, 3)
top10$PValue <- signif(top10$PValue, 3)

kable(top10,
      caption = "Table 1: Top 10 Differentially Expressed Genes")
```

## Top genes upregulated in small intestine

```{r}
si_genes <- results[results$logFC > 0, ]

si_sorted <- si_genes[order(si_genes$FDR), ]

top10_si <- head(si_sorted, 10)

top10_si$logFC <- round(top10_si$logFC, 2)
top10_si$logCPM <- round(top10_si$logCPM, 2)
top10_si$FDR <- signif(top10_si$FDR, 3)
top10_si$PValue <- signif(top10_si$PValue, 3)

kable(top10_si,
      caption = "Table 2: Top 10 Genes Upregulated in Small Intestine")
```

## Save full results table

```{r}
write.table(results,
            file = "small_intestine_vs_liver_edgeR_results.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)
```
