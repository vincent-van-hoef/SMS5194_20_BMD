---
jupytext:
    formats: md:myst
    text_representation:
      extension: .md
      format_name: myst
kernelspec:
      display_name: R
      language: R
      name: ir
---

# Analysis

Here you can find the results for project #5194 and its extension.

## Previous Work

The previous work for this project is included here as well.

```{code-cell} R
:tags: [hide-input]

setwd("/home/rstudio/scripts")

suppressWarnings(suppressMessages(suppressPackageStartupMessages({
library("kableExtra")
library("openxlsx")
library("ggplot2")
library("PCAtools")
library("factoextra")
library("FactoMineR")
library("cowplot")
library("ggpubr")
library("limma")
library("openxlsx")
library("matrixStats")
library("pcaExplorer")
library("dplyr")
})))
```


Load the user-supplied data.

```{code-cell} R
# Load and clean data...
MQ <- readRDS("../data/external/MQ.RData")
MNQ <- readRDS("../data/external/MNQ.RData")
MNQC <- readRDS("../data/external/MNQC.RData")
annotation <- readRDS("../data/external/annotation.RDS")
ann <- annotation[match(rownames(MNQC), annotation$Name),c(1:4, 12:19, 22:ncol(annotation))]
QCT <- read.csv("../data/external/SuperG_pQCT_PE.csv", sep = ";", dec=".", strip.white = TRUE)
phenoData <- readRDS("../data/external/phenoData.RData")
pheno	<- phenoData %>%
    as.data.frame() %>%
    select(-c(Basename, filenames, xMed, yMed, predictedSex)) %>%
    mutate(Fx = ifelse(Fx == 0, "noFx", "Fx")) %>%
    mutate(Batch = ifelse(Batch == 1, "Batch1", "Batch2")) %>%
    mutate(slide_array = paste(Slide, Array, sep="_")) %>%
    mutate(group = paste(BMD, Fx, sep="_")) %>%
    mutate(sample_name = gsub("blood.*", "", sample_name)) %>%
    mutate(sample_names= paste0("MKB", sample_name)) %>%
    left_join(select(QCT, 
                    id, 
                    TRAB_A_RAD_4, 
                    TRAB_DEN_RAD_4, 
                    TOT_A_RAD_4, 
                    TOT_DEN_RAD_4, 
                    CRT_A_RAD_4, 
                    CRT_DEN_RAD_4, 
                    CRT_THK_C_RAD_4), by = c("sample_names" = "id")) %>%
    select(-c(Slide, Array)) %>%
    mutate_at(.vars = vars("Chip.number", "Batch", "BMD", "Fx", "group"), factor) %>%
    tibble::column_to_rownames(var = "sample_name")

colnames(MQ)	<- rownames(pheno[match(colnames(MQ), pheno$slide_array),]) 
colnames(MNQ)	<- rownames(pheno[match(colnames(MNQ), pheno$slide_array),]) 
```


```{code-cell} R
m_all <- MNQ[order(rowVars(MNQ), decreasing = TRUE), ][,]
m_all <- m_all[,rownames(pheno)]
p_all <- prcomp(t(m_all), scale. = FALSE)
fviz_pca_ind(p_all,
         axes = c(1,2),
         geom.ind = c("point", "text"),
         pointshape = 21,
         fill.ind = pheno$group,
         palette = "jco",
         label = "",
         title = "All probes")
```



Digging deeper, we can list the p values of the correlation of the different types of metadata with the principal components (see Table \@ref(tab:pca-table)) to detect structure in the data. A Kruskal-Wallis (for categorical metadata) or a Spearman correlation (for continuous metadata) is performed on the first 10 principal components. Significant structure in the data would show up as significant p values.

```{code-cell} R
corr_all <- correlatePCs(p_all, coldata=pheno[, c("group", "Fx", "BMD", "Batch", "TRAB_A_RAD_4", "TRAB_DEN_RAD_4", "TOT_A_RAD_4", "TOT_DEN_RAD_4", "CRT_A_RAD_4", "CRT_DEN_RAD_4", "CRT_THK_C_RAD_4")], pcs = 1:10)
corr_all %>%
   as.data.frame(corr_all, row.names = rownames(corr_all)) %>%
   tibble::rownames_to_column() %>%
   mutate(rowname = gsub("_", "-", rowname)) %>%
   rename_all(function(x) gsub("_", " ", x)) %>%
   mutate_if(is.numeric, ~ round(., 4)) %>%
   mutate_if(is.numeric, ~cell_spec(., color = ifelse(. < 0.05, "red", "black"))) %>%
   tibble::column_to_rownames() %>%
   kable(escape = FALSE, row.names = TRUE, caption = "P values of correlation test of metadata with PCs. Kruskal-Willis test for categorical data and Spearman correlation for continuous data. PCA calculated using all probes.") %>%
   kable_styling(full_width = FALSE, latex_options = c("scale_down", "hold_position"))

```

The largest structure seems to be the Batch effect. This can be visualized using in particular PC4 vs PC5 in Figure \@ref(fig:pca-batch).

```{code-cell} R
fviz_pca_ind(p_all,
         axes = c(4,5),
         geom.ind = c("point", "text"),
         pointshape = 21,
         fill.ind = pheno$Batch,
         palette = "jco",
         label = "",
         title = "All probes")
```

Do we see more structure if we restrict the data to the 10 000 most variable probes? See Figure \@ref(fig:pcafig-var) for the PCA visualization using 10 000 most variable probes. Again, obvious structure in the data according to the biological conditions seems absent.

```{code-cell} R
m_var	<- m_all[1:10000,]
p_var <- prcomp(t(m_var), scale. = FALSE)
fviz_pca_ind(p_var,
         axes = c(1,2),
         geom.ind = c("point", "text"),
         pointshape = 21,
         fill.ind = pheno$group,
         palette = "jco",
         label = "",
         title = "10 000 most variable probes")
```

Again, we can list the p values of the correlation of the metadata with the principal components calculated using the most variable probes (see Table \@ref(tab:pca-var-table)). A Kruskal-Wallis (for categorical metadata) or a Spearman correlation (for continuous metadata) is performed on the first 10 principal components. Significant structure in the data would show up as significant p values.

```{code-cell} R
corr_var <- correlatePCs(p_var, coldata=pheno[, c("group", "Fx", "BMD", "Batch", "TRAB_A_RAD_4", "TRAB_DEN_RAD_4", "TOT_A_RAD_4", "TOT_DEN_RAD_4", "CRT_A_RAD_4", "CRT_DEN_RAD_4", "CRT_THK_C_RAD_4")], pcs = 1:10)
corr_var %>%
   as.data.frame(corr_var, row.names = rownames(corr_var)) %>%
   tibble::rownames_to_column() %>%
   mutate(rowname = gsub("_", "-", rowname)) %>%
   rename_all(function(x) gsub("_", " ", x)) %>%
   mutate_if(is.numeric, ~ round(., 4)) %>%
   mutate_if(is.numeric, ~cell_spec(., color = ifelse(. < 0.05, "red", "black"))) %>%
   tibble::column_to_rownames() %>%
   kable(escape = FALSE, row.names = TRUE, caption = "P values of correlation test of metadata with PCs. Kruskal-Willis test for categorical data and Spearman correlation for continuous data. PCA calculated using the 10 000 most variable probes.") %>%
   kable_styling(full_width = FALSE, latex_options = c("scale_down", "hold_position"))
```

The batch effect is even a little bit more clear when looking at the PC2 vs PC3 of this most variable subset of the data.

```{code-cell} R
fviz_pca_ind(p_var,
         axes = c(2,3),
         geom.ind = c("point", "text"),
         pointshape = 21,
         fill.ind = pheno$Batch,
         palette = "jco",
         label = "",
         title = "10 000 most variable probes")
```

CONCLUSION: There is not so much structure in the data except for the batch effect. There seems to be some correlation with the continuous phenotypical covariates. How does this reflect on the differential methylation?

## Differential methylation BMD and Fractures

Earlier analysis was done on Combat adjusted data and did not result in significantly differentially methylated probes for the Fracture or BMD comparisons. The question now is whether including the batch effect in the linear model might yield more results.

### Probe Level

To assess the different methylation levels, we used the limma R package. Covariates included in each model are the batch effect and the cell type composition as calculated by the group.

After multiple correction no single probes pass the significance threshold. A quick literature search reveals that on the single probe level several bone/osteoporosis articles failed to show significant methylation levels.

Full results for these differential methylation comparisons are collected in the "Probes" folder.

### Region Level

Differential methylation is often evaluated at a regional level - so not single changed CpG sites but multiple close sites that are methylated similarly.

[mCGSEA (methylated CpGs Set Enrichment Analysis)](https://pubmed.ncbi.nlm.nih.gov/30753302/) is a GSEA-based differential methylation analysis where gene sets are defined as sets of CpG sites in predefined regions. This new tool is capable to detect subtle but consistent methylation differences in predefined genomic regions from 450 K and EPIC microarrays data. The predefined regions are promoter regions, genes and CpG islands.

**mCSEA is based on Gene Set Enrichment analysis (GSEA), a popular methodology for functional analysis that was specifically designed to avoid some drawbacks in the field of gene expression. GSEA is able to detect significant gene sets that exhibit strong cross-correlation when differential expression of individual genes is modest from the statistical point of view. GSEA uses a given statistical metric to rank all genes of a genome and applies a weighted Kolmogorov–Smirnov (KS) statistic to calculate an Enrichment Score (ES). Basically, ES for each set is calculated running through the entire ranked list increasing the score when a gene in the set is encountered and decreasing the score when the gene encountered is not in the analyzed set. ES of this set is the maximum difference from 0.**

The statistic that is used in this package to rank the probes is the t-statistic from the limma test of the individual probes. This statistic is somewhat similar to the log2 fold change normalized to the standard error of the probe.

The ES and NES (ES normalized for set size) of a set are calculated by looking at where the statistics of probes belonging to a certain set can be found in the ranked probe list. A high (N)ES indicates these probes are found high up in the list of "BMD high vs low" or "fractures vs no fractures". In other words, a high (N)ES value means that for the probes in this set there is - on average - a shift towards a higher methylation in high BMD or fracture versus no fracture.

The significance of each (N)ES is calculated permuting the sets and recomputing ES, getting a null distribution for the ES. A multiple comparison correction is also performed on the p values.

Results for these comparisons are collected in the "Regions" folder. Each excel result file has three tabs for the three predefined regions:

* Genes: probes belonging to the body of a gene
* Promoters: probes belonging to the promoter of a gene
* CpG Islands: probes belonging to CpG Islands

More info on method and results can be found in the linked paper at the top of this section.

```{code-cell} R

mains	<- c("BMD", "Fx")

for(main in mains){
mod <- model.matrix(as.formula(paste0("~ 0 + ", main, " + Batch + CD8T + CD4T + NK + Bcell + Mono + Neu")), data = pheno)
fit <- lmFit(MNQ, mod)
if(main == "Fx"){
   contMatrix <- makeContrasts(contrasts = "FxFx-FxnoFx", levels = mod)
} else {
   contMatrix <- makeContrasts(contrasts = "BMDHigh-BMDLow", levels = mod)
}
fit_contrast	<- contrasts.fit(fit, contMatrix)
fit_res	<- eBayes(fit_contrast)

#write.xlsx(topTable(fit_res, n=Inf), paste0("Results/Probes/", main, ".xlsx"), row.names=TRUE)

myRank 	<- rankProbes(MNQ,
                      pheno, 
                      explanatory = main, 
                      covariates = c("Batch", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                      continuous = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), 
                      typeInput = "M")
myTest	<-  mCSEATest(myRank, 
                      MNQ, 
                      pheno = pheno, 
                      column = main, 
                      regionsTypes = c("promoters", "genes", "CGI"), 
                      platform = "EPIC")

#write.xlsx(myTest[c("genes", "promoters", "CGI")], paste0("Results/Regions/mCSEA_", main, ".xlsx"), row.names = TRUE) 
}
```

## Differential methylation QCT

### Probe Level

To assess the different methylation levels for the QCT measurements, we again used the limma R package. Covariates included in each model are again the batch effect and the cell type composition as calculated by the group.

After multiple correction no single probes pass the significance threshold.

As before, results for these comparisons are collected in the "Probes" folder.

### Region Level

As before, probes are ranked according to the t-statistic of their association with the QCT value and used as input for the mCSEA analysis. A higher (N)ES value means that for the probes in this set there is - on average - a shift towards a higher methylation in higher values of the QCT measurement.

Results for these comparisons are also collected in the "Regions" folder. Each excel has three tabs:

* Genes: probes belonging to the body of a gene
* Promoters: probes belonging to the promoter of a gene
* CpG Islands: probes belonging to CpG Islands

```{code-cell} R

QCTs	<- c("TRAB_A_RAD_4",
            "TRAB_DEN_RAD_4",
            "TOT_A_RAD_4",
            "TOT_DEN_RAD_4",
            "CRT_A_RAD_4",
            "CRT_DEN_RAD_4",
            "CRT_THK_C_RAD_4")

for(qct in QCTs){
mod <- model.matrix(as.formula(paste0("~ ", qct, " + Batch + CD8T + CD4T + NK + Bcell + Mono + Neu")), data = pheno)
M_cont	<- MNQ[, rownames(mod)]
fit <- lmFit(M_cont, mod)
fit_contrast 	<- contrasts.fit(fit, coefficients = 2)
fit_res	<- eBayes(fit_contrast)

#write.xlsx(topTable(fit_res, n=Inf), paste0("Results/Probes/", qct, ".xlsx"), row.names=TRUE)

myRank 	<- setNames(as.numeric(fit_res$t), rownames(fit_res$t))
myTest	<-  mCSEATest(myRank, 
                      MNQ, 
                      pheno = pheno, 
                      column = qct, 
                      regionsTypes = c("promoters", "genes", "CGI"), 
                      platform = "EPIC")
#write.xlsx(myTest[c("genes", "promoters", "CGI")], paste0("Results/Regions/mCSEA_", qct, ".xlsx"), row.names = TRUE) 
}
```

## QCT Data

### PCA

To look for further structure in the QCT data, a PCA on this data was performed as well (see Figure \@ref(fig:qct-pca)). On visual inspection, it seems there might be a slight difference between Low_NoFx and the other groups.

```{code-cell} R
# Perform a pca pn the clinical data only
pheno_clin 	<- pheno %>%
  as.data.frame() %>%
  dplyr::select(ends_with("4")) %>%
  drop_na() %>%
  as.matrix()
p_clin <- prcomp(pheno_clin, scale. = TRUE)
fviz_pca_ind(p_clin,
            axes = c(1,2),
            #geom.var = c("arrow", "text"),
            pointshape = 21,
            geom.ind = "point",
            fill.ind = pheno[rownames(pheno_clin), "group"],
            legend.title = "Group",
            #repel = TRUE,
            title = "")
```

### Comparative statistics

Simple comparative statistics per group for the different QCT measurements are visualized below and give a bit more info on the data structure.

```{code-cell} R
dat 	<- pheno %>%
  pivot_longer(cols = ends_with("RAD_4"),
  names_to = "metric",
  values_to = "value",
  values_drop_na = TRUE)
ggboxplot(dat, x = "BMD", y = "value", color = "BMD", add = "jitter") +
                facet_wrap( ~ metric, scales = "free") +
                stat_compare_means(method = "t.test") +
                scale_y_continuous(expand = expansion(mult = 0.1)) 
```

```{code-cell} R
dat 	<- pheno %>%
  pivot_longer(cols = ends_with("RAD_4"),
  names_to = "metric",
  values_to = "value",
  values_drop_na = TRUE)
ggboxplot(dat, x = "Fx", y = "value", color = "Fx", add = "jitter") +
                facet_wrap( ~ metric, scales = "free") +
                stat_compare_means(method = "t.test") +
                scale_y_continuous(expand = expansion(mult = 0.1))
```

```{code-cell} R
dat 	<- pheno %>%
  pivot_longer(cols = ends_with("RAD_4"),
  names_to = "metric",
  values_to = "value",
  values_drop_na = TRUE)
ggboxplot(dat, x = "group", y = "value", color = "group", add = "jitter") +
                facet_wrap( ~ metric, scales = "free") +
                stat_compare_means(method = "anova") +
                scale_y_continuous(expand = expansion(mult = 0.1))
```


## Extension