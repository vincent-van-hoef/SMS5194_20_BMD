suppressWarnings(suppressMessages(suppressPackageStartupMessages({
library("kableExtra")
library("knitr")
library("openxlsx")
library("ggplot2")
library("PCAtools")
library("factoextra")
library("FactoMineR")
library("cowplot")
library("ggpubr")
library("limma")
library("openxlsx")
library("dplyr")
})))

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

corr_all <- correlatePCs(p_all, coldata=pheno[, c("group", "Fx", "BMD", "Batch", "TRAB_A_RAD_4", "TRAB_DEN_RAD_4", "TOT_A_RAD_4", "TOT_DEN_RAD_4", "CRT_A_RAD_4", "CRT_DEN_RAD_4", "CRT_THK_C_RAD_4")], pcs = 1:10)
corr_all %>%
   as.data.frame(corr_all, row.names = rownames(corr_all)) %>%
   rownames_to_column() %>%
   mutate(rowname = gsub("_", "-", rowname)) %>%
   rename_all(function(x) gsub("_", " ", x)) %>%
   mutate_if(is.numeric, ~ round(., 4)) %>%
   mutate_if(is.numeric, ~cell_spec(., color = ifelse(. < 0.05, "red", "black"))) %>%
   column_to_rownames() %>%
   kable(escape = FALSE, row.names = TRUE, caption = "P values of correlation test of metadata with PCs. Kruskal-Willis test for categorical data and Spearman correlation for continuous data. PCA calculated using all probes.") %>%
   kable_styling(full_width = FALSE, latex_options = c("scale_down", "hold_position"))


fviz_pca_ind(p_all,
         axes = c(4,5),
         geom.ind = c("point", "text"),
         pointshape = 21,
         fill.ind = pheno$Batch,
         palette = "jco",
         label = "",
         title = "All probes")

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

corr_var <- correlatePCs(p_var, coldata=pheno[, c("group", "Fx", "BMD", "Batch", "TRAB_A_RAD_4", "TRAB_DEN_RAD_4", "TOT_A_RAD_4", "TOT_DEN_RAD_4", "CRT_A_RAD_4", "CRT_DEN_RAD_4", "CRT_THK_C_RAD_4")], pcs = 1:10)
corr_var %>%
   as.data.frame(corr_var, row.names = rownames(corr_var)) %>%
   rownames_to_column() %>%
   mutate(rowname = gsub("_", "-", rowname)) %>%
   rename_all(function(x) gsub("_", " ", x)) %>%
   mutate_if(is.numeric, ~ round(., 4)) %>%
   mutate_if(is.numeric, ~cell_spec(., color = ifelse(. < 0.05, "red", "black"))) %>%
   column_to_rownames() %>%
   kable(escape = FALSE, row.names = TRUE, caption = "P values of correlation test of metadata with PCs. Kruskal-Willis test for categorical data and Spearman correlation for continuous data. PCA calculated using the 10 000 most variable probes.") %>%
   kable_styling(full_width = FALSE, latex_options = c("scale_down", "hold_position"))

fviz_pca_ind(p_var,
         axes = c(2,3),
         geom.ind = c("point", "text"),
         pointshape = 21,
         fill.ind = pheno$Batch,
         palette = "jco",
         label = "",
         title = "10 000 most variable probes")


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

write.xlsx(topTable(fit_res, n=Inf), paste0("Results/Probes/", main, ".xlsx"), row.names=TRUE)

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

write.xlsx(myTest[c("genes", "promoters", "CGI")], paste0("Results/Regions/mCSEA_", main, ".xlsx"), row.names = TRUE) 
}


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

write.xlsx(topTable(fit_res, n=Inf), paste0("Results/Probes/", qct, ".xlsx"), row.names=TRUE)

myRank 	<- setNames(as.numeric(fit_res$t), rownames(fit_res$t))
myTest	<-  mCSEATest(myRank, 
                      MNQ, 
                      pheno = pheno, 
                      column = qct, 
                      regionsTypes = c("promoters", "genes", "CGI"), 
                      platform = "EPIC")
write.xlsx(myTest[c("genes", "promoters", "CGI")], paste0("Results/Regions/mCSEA_", qct, ".xlsx"), row.names = TRUE) 
}

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

dat 	<- pheno %>%
  pivot_longer(cols = ends_with("RAD_4"),
  names_to = "metric",
  values_to = "value",
  values_drop_na = TRUE)
ggboxplot(dat, x = "BMD", y = "value", color = "BMD", add = "jitter") +
                facet_wrap( ~ metric, scales = "free") +
                stat_compare_means(method = "t.test") +
                scale_y_continuous(expand = expansion(mult = 0.1)) 

dat 	<- pheno %>%
  pivot_longer(cols = ends_with("RAD_4"),
  names_to = "metric",
  values_to = "value",
  values_drop_na = TRUE)
ggboxplot(dat, x = "Fx", y = "value", color = "Fx", add = "jitter") +
                facet_wrap( ~ metric, scales = "free") +
                stat_compare_means(method = "t.test") +
                scale_y_continuous(expand = expansion(mult = 0.1))

dat 	<- pheno %>%
  pivot_longer(cols = ends_with("RAD_4"),
  names_to = "metric",
  values_to = "value",
  values_drop_na = TRUE)
ggboxplot(dat, x = "group", y = "value", color = "group", add = "jitter") +
                facet_wrap( ~ metric, scales = "free") +
                stat_compare_means(method = "anova") +
                scale_y_continuous(expand = expansion(mult = 0.1))
