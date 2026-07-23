library(stringr)
library(tidyverse)
library(celda)
library(singleCellTK)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)
library(scater)
library(nlme)
library(ggtext)
library(GSVA)
library(biomaRt)
library(limma)
library(hypeR)
library(emmeans)
library(SingleCellExperiment)
library(ggpmisc)
library(caret)
library(e1071)
library(pROC)

############################## Set seed
# set.seed(123)

################################################################################
######### Read in SCE
sce <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/Fig4-6_Paper_Scripts_Results/Misc_Scripts/EpithelialBrushBiopsyCells_SCE_Final.rds")
sce$epithelial_CellType <- gsub("^Basal Cells", "LG Basal Cells", sce$epithelial_CellType)
sce_orig <- sce

celltype_colors <- readRDS("../celltype_colors.rds")
smoking_colors <- readRDS("../smoking_colors.rds")
sampletype_cols <- readRDS("../sampletype_cols.rds")
names(sampletype_cols)[1] <- "Biopsy"
histology_colors <- readRDS("../histology_colors.rds")

##############################################################################################################
##################################### Figure 6A

biopsy_epi <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/NonimmuneBiopsyCell_SCE_Final.rds")

#### Get column names for biopsy cells
shared_cells <- intersect(colnames(biopsy_epi), colnames(sce))

#### Map cell names to column indices
sce_idx <- match(shared_cells, colnames(sce))
biopsy_idx <- match(shared_cells, colnames(biopsy_epi))

#### Initialize empty metadata columns in sce
sce$biopsy_celda_cluster <- NA
sce$biopsy_epithelial_CellType <- NA
sce$biopsy_CellType_Cluster <- NA

#### Transfer metadata using numeric indices
sce$biopsy_celda_cluster[sce_idx] <- biopsy_epi$celda_cell_cluster[biopsy_idx]
sce$biopsy_epithelial_CellType[sce_idx] <- biopsy_epi$CellType[biopsy_idx]
sce$biopsy_CellType_Cluster[sce_idx] <- biopsy_epi$CellType_Cluster[biopsy_idx]

#### Prep UMAP
altsce <- altExp(sce)
umap_celda = altsce@int_colData@listData[["reducedDims"]]@listData[["celda_UMAP"]] %>% 
  as.data.frame() %>% cbind(biopsy_celda_cluster = sce$biopsy_celda_cluster, 
                            celda_cluster = sce$epithelial_celda_cluster)

umap_celda <- umap_celda %>%
  mutate(legend_category = case_when(
    biopsy_celda_cluster %in% 1:6 ~ "HG Basal Cells",
    biopsy_celda_cluster %in% 7:11 ~ "LG Basal Cells",
    celda_cluster == 16 ~ "HG Basal Brush Cells",
    TRUE ~ NA_character_  # Assign NA to other points if not part of these categories
  ))

############ Create & plot UMAP
pdf("./Figures/Fig6A_HGB_vs_LGB_BiopsyCells.pdf", width = 8, height = 7)
ggplot() +
  # Plot "other" points first in grey90, so they appear underneath
  geom_point(data = umap_celda %>% filter(is.na(legend_category)),
             aes(x = celda_UMAP1, y = celda_UMAP2), color = "grey90", size = 2) +
  
  # Plot "HG Basal Cells" next
  geom_point(data = umap_celda %>% filter(legend_category == "HG Basal Cells"),
             aes(x = celda_UMAP1, y = celda_UMAP2, color = legend_category), size = 2) +
  
  # Plot "Basal Cells" next
  geom_point(data = umap_celda %>% filter(legend_category == "LG Basal Cells"),
             aes(x = celda_UMAP1, y = celda_UMAP2, color = legend_category), size = 2) +
  
  # # Plot "HG Basal Brush Cells" last to ensure it appears on top
  # geom_point(data = umap_celda %>% filter(legend_category == "HG Basal Brush Cells"),
  #            aes(x = celda_UMAP1, y = celda_UMAP2, color = legend_category), size = 2) +
  
  # Define color scheme for the legend categories
  scale_color_manual(
    values = c(
      "HG Basal Cells" = "red",
      "LG Basal Cells" = "blue"
      # "HG Basal Brush Cells" = "#FFC20A"
    ),
    na.value = "grey90"  # Ensures other points are still plotted in grey90
  ) +
  
  # Customize theme and set legend title
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(color = "Basal Cell Subtype")  # Set legend title
dev.off()

##############################################################################################################
##################################### Figure 6B & 6C

sce <- sce_orig 

sce$Subject_TimePoint_SampleType <- paste0(sce$Subject, "_", sce$TimePoint, "_", sce$SampleType)

hgb <- sce[,which(sce$epithelial_CellType == "HG Basal Cells")]
table(hgb$SampleType, hgb$Subject_TimePoint_SampleType)

subj152 <- hgb[,which(hgb$Subject == "152")]
table(subj152$Subject_TimePoint_SampleType, subj152$Anatomic_Site)

################### Plot Subject 152 at T6 by Sample Type
altsce <- altExp(sce)
umap_celda = altsce@int_colData@listData[["reducedDims"]]@listData[["celda_UMAP"]] %>% 
  as.data.frame() %>% cbind(sampletype = sce$Subject_TimePoint_SampleType)

legend_title = "Sample Types at T6"

pdf("./Figures/Fig6B_Subj152_SampleType_Final.pdf")
ggplot(umap_celda, aes(x=celda_UMAP1, y=celda_UMAP2, color=sampletype)) + geom_point() + 
  scale_color_manual(legend_title, values=c("152_T6_Biopsy" = "#0C7BDC",
                                            # "PCGA02_20152_3100046_1" = "red",
                                            # "PCGA02_20152_3100054_1" = "#0C7BDC", 
                                            "152_T6_BronchialBrush" = "#FFC20A"),
                     na.value = "grey90",
                     labels = c("Biopsy", "BronchialBrush")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # guides(fill=guide_legend(title="New Legend Title")) +
  xlab("UMAP1") + ylab("UMAP2")
dev.off()

################### Plot Subject 152 at T6 by Sample
altsce <- altExp(sce)
umap_celda = altsce@int_colData@listData[["reducedDims"]]@listData[["celda_UMAP"]] %>% 
  as.data.frame() %>% cbind(sampletype = sce$PCGA02_CellID_Template)

legend_title = "Samples Taken at T6"

pdf("./Figures/Fig6C_Subj152_Samples_Final.pdf")
ggplot(umap_celda, aes(x=celda_UMAP1, y=celda_UMAP2, color=sampletype)) + geom_point() + 
  scale_color_manual(legend_title, values=c("PCGA02_20152_3100053_1" = "red",
                                            # "PCGA02_20152_3100046_1" = "red",
                                            "PCGA02_20152_3100054_1" = "blue", 
                                            "PCGA02_20152_3600055_1" = "#FFC20A"),
                     na.value = "grey90",
                     labels = c("LUSC Biopsy at RB4/5", "CIS Biopsy at LMB", "Bronchial Brush at LB6")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # guides(fill=guide_legend(title="New Legend Title")) +
  xlab("UMAP1") + ylab("UMAP2")
dev.off()

##############################################################################################################
##################################### Figure 6D

sce <- sce_orig 

########## Identify modules differentially expressed between clusters of the same cell type
f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")), 
                                    colData = colData(sce))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")

sce_modular <- sce_modular[,which(sce_modular$SampleType == "BronchialBrush")]

cellprobs <- data.frame(t(assay(sce_modular,"LogNormalize")))
cellprobs$Cluster <- sce_modular$epithelial_celda_cluster ## 
cellprobs$CellType <- sce_modular$epithelial_CellType
cellprobs$Sample <- sce_modular$Sample
cellprobs$Subject <- sce_modular$Subject
cellprobs$Smoking_Status <- sce_modular$Smoking_Status
cellprobs$SampleType <- sce_modular$SampleType
cellprobs$Site <- sce_modular$PCGA02_site

epi_clust_de <- function(expr = cellprobs, group_var = "CellType", clusters) {
  expr_subset <- expr
  each_model_res <- lapply(clusters, function(x) {
    expr_subset[,"EpiCluster"] <- ifelse(expr_subset[,group_var] == x, "B", "A")
    model_p <- data.frame(Cluster = x, Module = names(expr_subset)[1:70], P = rep(0, 70), Beta = rep(0, 70), Group_Mean = rep(0, 70), NotGroup_Mean = rep(0, 70), row.names = names(expr_subset)[1:70]) # Smoking_P = rep(0, 55), Smoking_Beta = rep(0, 55), 
    for(i in rownames(model_p)) {
      model_p[i,"Group_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] == x),i], na.rm = TRUE)
      model_p[i,"NotGroup_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] != x),i], na.rm = TRUE)
      
      # model_construct <- as.formula(paste0(i," ~ EpiCluster"))
      if(group_var == "SampleType") {
        model_construct <- as.formula(paste0(i," ~ EpiCluster"))
      } else {
        model_construct<- as.formula(paste0(i," ~ EpiCluster + Smoking_Status + Site"))
      }
      model <- lme(model_construct,
                   random = ~ 1|Sample,
                   data = expr_subset,
                   control = lmeControl(opt = "optim"),
                   # control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
                   na.action = na.omit)
      tTab <- summary(model)$tTable
      rownames(tTab)[2] <- "EpiCluster"
      model_p[i,"Beta"] <- tTab["EpiCluster","Value"]
      model_p[i,"P"] <- tTab["EpiCluster","p-value"]
    }
    model_p$Q <- p.adjust(model_p$P, n = nrow(model_p), method = "fdr")
    model_p$LogQ <- -log10(model_p$Q)
    model_p$Significant <- ((model_p$Q < 0.05) & (abs(model_p$Beta) > 0.5))
    model_p$SignificantPositive <- ((model_p$Q < 0.05) & (model_p$Beta > 0.5))
    model_p <- model_p %>% arrange(desc(SignificantPositive), Q)
  })
  names(each_model_res) <- clusters
  each_model_res
}

################### Find modules upregulated in HGB vs LGB bronchial brush cells

# basal_cellprobs <- cellprobs[which(cellprobs$CellType %in% c("LG Basal Cells", "HG Basal Cells")),]
# diff_modules_basal <- epi_clust_de(expr = basal_cellprobs, clusters = unique(basal_cellprobs$CellType))
# diff_modules_res_basal <- Reduce(diff_modules_basal, f = rbind)
# diff_modules_res_significant_basal <- diff_modules_res_basal[which(diff_modules_res_basal$SignificantPositive),]

# saveRDS(diff_modules_res_significant_basal, "./Objects_Results/Fig6D_DEGModules_Significant_Final.rds")

diff_modules_res_significant_basal <- readRDS("./Objects_Results/Fig6D_DEGModules_Significant_Final.rds")

############ Basal Heatmap
modules <- unique(diff_modules_res_significant_basal$Module)

########## The factorizeMatrix function generates a "factorized" matrix between features/feature modules and samples/cell clusters/cells.
factoMat <- factorizeMatrix(sce, useAssay = "decontXcounts")
cellModuleMatrix <- factoMat$counts$cell

########## Standard normalization of counts
pop <- SingleCellExperiment(assays = SimpleList("factorizeCounts" = cellModuleMatrix), colData = colData(sce))
pop <- logNormCounts(pop, assay.type = "factorizeCounts", name = "logcounts")

pop <- pop[,which(pop$epithelial_CellType %in% c("LG Basal Cells", "HG Basal Cells"))]
pop <- pop[,which(pop$SampleType %in% c("BronchialBrush"))]

########## LogFC threshold
logFC <- as.data.frame(matrix(NA, nrow = length(modules)), ncol = 2)
logFC[,1] <- modules

for (i in 1:length(modules)) {
  mod_nas <- mean(logcounts(pop)[modules[i],which(pop$epithelial_CellType == "HG Basal Cells")])
  mod_other <- mean(logcounts(pop)[modules[i],which(!pop$epithelial_CellType == "HG Basal Cells")])
  logFC[i,2] <- mod_nas - mod_other
}
colnames(logFC) <- c("Module", "LogFC")
modules <- logFC$Module[which(abs(logFC$LogFC) > 1.5)]

diff_modules_res_significant_basal <- diff_modules_res_significant_basal[which(diff_modules_res_significant_basal$Module %in% modules),]
# saveRDS(diff_modules_res_significant_basal, "./Objects_Results/Fig6D_DEGModules_Significant_LogFC_Final.rds")

cellModuleMatrix_norm <- as.data.frame(assay(pop, "logcounts"))

########## Subset factorized proportions matrix to just the modules of interest
popNorm2 <- cellModuleMatrix_norm[which(rownames(cellModuleMatrix_norm) %in% modules),]

########## Heatmap of individual celda module scores
popNorm3 <- t(scale(t(popNorm2)))

########## Reorder rows and columns & rename modules
popNorm3 <- popNorm3[c("L2", "L9", "L10", "L11", "L12", "L13", "L14", "L15", "L21", "L33", "L63", "L64",
                       "L59", "L62", "L66"),]

rownames(popNorm3)[which(rownames(popNorm3) == "L2")] <- "L2 (CD74)"
rownames(popNorm3)[which(rownames(popNorm3) == "L9")] <- "L9 (SCGB1A1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L10")] <- "L10 (SCGB3A1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L11")] <- "L11 (BPIFA1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L12")] <- "L12 (BPIFB1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L13")] <- "L13 (SLPI, C3)"
rownames(popNorm3)[which(rownames(popNorm3) == "L14")] <- "L14 (MSMB)"
rownames(popNorm3)[which(rownames(popNorm3) == "L15")] <- "L15 (MUC5AC)"
rownames(popNorm3)[which(rownames(popNorm3) == "L21")] <- "L21 (WFDC2, PIGR)"
rownames(popNorm3)[which(rownames(popNorm3) == "L33")] <- "L33 (EPAS1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L63")] <- "L63 (KRT15)"
rownames(popNorm3)[which(rownames(popNorm3) == "L64")] <- "L64 (FOS, JUN)"
rownames(popNorm3)[which(rownames(popNorm3) == "L59")] <- "L59 (S100A9, SPRR2A)"
rownames(popNorm3)[which(rownames(popNorm3) == "L62")] <- "L62 (KRT6A)"
rownames(popNorm3)[which(rownames(popNorm3) == "L66")] <- "L66 (GSTM2)"

########## Prep heatmap annotations
pop$epithelial_CellType <- factor(pop$epithelial_CellType, levels = c("LG Basal Cells", "HG Basal Cells"))

column_ha <- HeatmapAnnotation(Smoking_Status = pop$Smoking_Status,
                               Histology = pop$Worst_Histology_at_TimePoint_Label,
                               Cell_Type = pop$epithelial_CellType,
                               # Celda_Cluster = pop$epithelial_celda_cluster,
                               col = list(Smoking_Status = smoking_colors,
                                          # SampleType = sampletype_cols,
                                          Histology = histology_colors,
                                          Cell_Type = celltype_colors),
                               # Celda_Cluster = celda_cluster_cols),
                               show_annotation_name = TRUE,
                               show_legend = TRUE,
                               annotation_name_rot = 0)

popNorm3[which(popNorm3 > 2)] <- 2
popNorm3[which(popNorm3 < -2)] <- -2

########## Plot & save heatmap
pdf("./Figures/Fig6D_BasalModules_Heatmap_Final.pdf", width = 9, height = 7)
print(Heatmap(popNorm3, 
              top_annotation = column_ha, 
              show_column_names = FALSE, 
              column_split = pop$epithelial_CellType,
              show_heatmap_legend = TRUE,
              #column_title = NULL,
              # column_title_rot = 45,
              column_title_gp = gpar(fontsize = 10),
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()


##############################################################################################################
##################################### Figure 6E

##################### Read in data
t <- read.table("../Figure5/072324_Biopsies+Brushes_noIonocytes_decontXcounts_zCellType_proteinCodingGenes_module_features_K40_L70.csv", sep = "\t", header = TRUE)
t2 <- as.list(t)
t3 <- lapply(t2, function(z){ z[!is.na(z) & z != ""]})
names(t3) <- gsub("V", "L", names(t3))
# basal_modules <- t3[which(names(t3) %in% c("L33", "L62", "L70"))]
basal_modules <- t3[which(names(t3) %in% c("L2", "L10", "L12", "L13", "L14", "L15", "L21", "L33", "L64", "L59", "L62", "L66", "L70"))]

##################### PCGA1
resid_br_Disc <- readRDS("./Objects_Results/Fig6E_PCGA1_BronchialBrush_DiscoveryCohort_Residuals.rds")
resid_br_Val <- readRDS("./Objects_Results/Fig6E_PCGA1_BronchialBrush_ValidationCohort_Residuals.rds")

########## Running GSVA (already completed + saved)
### Convert ensembl IDs to gene symbols
# bulk_br <- readRDS("/restricted/projectnb/pcga/ANALYSIS_NEXTFLOW/PCGA_nextflow_br.rds")
# m <- match(rownames(resid_br_Disc),rownames(bulk_br))
# rownames(resid_br_Disc) <- rowData(bulk_br)$hgnc_symbol[m]
# m <- match(rownames(resid_br_Val),rownames(bulk_br))
# rownames(resid_br_Val) <- rowData(bulk_br)$hgnc_symbol[m]

### Run GSVA
# gsvaparam_disc <- gsvaParam(exprData = resid_br_Disc, geneSets = basal_modules, kcdf = "Gaussian", absRanking = FALSE)
# br.disc_gsva.scores <- gsva(param  = gsvaparam_disc, verbose=TRUE)
# gsvaparam_val <- gsvaParam(exprData = resid_br_Val, geneSets = basal_modules, kcdf = "Gaussian", absRanking = FALSE)
# br.val_gsva.scores <- gsva(param = gsvaparam_val, verbose=TRUE)

### Save RDS files
# saveRDS(br.disc_gsva.scores, "./Objects_Results/Fig6E_PCGA1_Discovery_GSVA_Scores.rds")
# saveRDS(br.val_gsva.scores, "./Objects_Results/Fig6E_PCGA1_Validation_GSVA_Scores.rds")

### Read in GSVA objects
br.disc_gsva.scores <- readRDS("./Objects_Results/Fig6E_PCGA1_Discovery_GSVA_Scores.rds")
br.val_gsva.scores <- readRDS("./Objects_Results/Fig6E_PCGA1_Validation_GSVA_Scores.rds")

### Read in bulk data: residuals
clinical_annot <- readRDS("/restricted/projectnb/pcga/ANALYSIS_NEXTFLOW/subtype_associations/samples_annotation_v3.rds")

### Subset annotations
annot_br_Disc <- clinical_annot[colnames(resid_br_Disc),]
annot_br_Val <- clinical_annot[colnames(resid_br_Val),]

### Bulk annotations
br.disc.df <- cbind(annot_br_Disc, t(br.disc_gsva.scores))
br.val.df <- cbind(annot_br_Val, t(br.val_gsva.scores))

### Combine objects
pcga1.br.df <- rbind(br.disc.df, br.val.df)

### Edit metadata
pcga1.br.df$Dysplasia_Grade <- gsub("MildD", "Mild Dysplasia", pcga1.br.df$Dysplasia_Grade)
pcga1.br.df$Dysplasia_Grade <- gsub("ModD", "Moderate Dysplasia", pcga1.br.df$Dysplasia_Grade)
pcga1.br.df$Dysplasia_Grade <- gsub("SevD", "Severe Dysplasia", pcga1.br.df$Dysplasia_Grade)

pcga1.br.df$Dysplasia_Grade <- factor(pcga1.br.df$Dysplasia_Grade,
                                      levels = c('Normal', 'Hyperplasia', 'Metaplasia', 'Mild Dysplasia', 'Moderate Dysplasia', 'Severe Dysplasia'))

pcga1.br.df$Dysplasia_Grade2 <- as.numeric(pcga1.br.df$Dysplasia_Grade)


############# Test if there is a positive linear trend

lmRes_L33_br <- lme(L33 ~ Dysplasia_Grade2 + L70 + Cohort + Genomic_Smoking_Status, random = ~ 1|Patient, data = pcga1.br.df, na.action = na.omit)
summary(lmRes_L33_br)
saveRDS(lmRes_L33_br, "./Objects_Results/Fig6E_PCGA1_ModuleL33_linearModelResults.rds")

lmRes_L62_br <- lme(L62 ~ Dysplasia_Grade2 + L70 + Cohort + Genomic_Smoking_Status, random = ~ 1|Patient, data = pcga1.br.df, na.action = na.omit)
summary(lmRes_L62_br)
saveRDS(lmRes_L62_br, "./Objects_Results/Fig6E_PCGA1_ModuleL62_linearModelResults.rds")

lmRes_L70_br <- lme(L70 ~ Dysplasia_Grade2 + Cohort + Genomic_Smoking_Status, random = ~ 1|Patient, data = pcga1.br.df, na.action = na.omit)
summary(lmRes_L70_br)
saveRDS(lmRes_L70_br, "./Objects_Results/Fig6E_PCGA1_ModuleL70_linearModelResults.rds")

############# Save results from all modules
# Initialize a list to store results
modules <- c("L2", "L10", "L12", "L13", "L14", "L15", "L21", "L33", "L64", "L59", "L62", "L66")
model_results <- list()

for (mod in modules) {
  # Build the formula dynamically
  formula_str <- paste0(mod, " ~ Dysplasia_Grade2 + L70 + Cohort + Genomic_Smoking_Status")
  formula_obj <- as.formula(formula_str)
  
  # Fit the model
  model <- lme(formula_obj, random = ~ 1 | Patient, data = pcga1.br.df, na.action = na.omit)
  
  # Store the summary in a named list
  model_results[[mod]] <- summary(model)
}

model_results[["L70"]] <- summary(lmRes_L70_br)

saveRDS(model_results, "./Objects_Results/Fig6E_PCGA1_All6DModules_linearModelResults.rds")

############################ Plotting

model <- lme(L33 ~ L70 + Cohort + Genomic_Smoking_Status, random = ~ 1|Patient, data = pcga1.br.df, na.action = na.omit)
residuals_L33 <- residuals(model)
pcga1.br.df$residuals_L33 <- residuals_L33

model <- lme(L62 ~ L70 + Cohort + Genomic_Smoking_Status, random = ~ 1|Patient, data = pcga1.br.df, na.action = na.omit)
residuals_L62 <- residuals(model)
pcga1.br.df$residuals_L62 <- residuals_L62

model <- lme(L70 ~ Cohort + Genomic_Smoking_Status, random = ~ 1|Patient, data = pcga1.br.df, na.action = na.omit)
residuals_L70 <- residuals(model)
pcga1.br.df$residuals_L70 <- residuals_L70

histology_colors <- readRDS("../histology_colors.rds")
names(histology_colors)[2] <- "Normal"

########## Module L33
pdf("./Figures/Fig6E_PCGA1_ModuleL33_Boxplot.pdf", height = 5, width = 4)
ggplot(pcga1.br.df, aes(x = Dysplasia_Grade, y = residuals_L33, fill = Dysplasia_Grade)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  scale_fill_manual(values = histology_colors) +
  theme_classic() +
  xlab("") + 
  ylab("L33 Residuals") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(fill = "none") +
  ggtitle("L33 (PCGA1)")
dev.off()

########## Module L62
pdf("./Figures/Fig6E_PCGA1_ModuleL62_Boxplot.pdf", height = 5, width = 4)
ggplot(pcga1.br.df, aes(x = Dysplasia_Grade, y = residuals_L62, fill = Dysplasia_Grade)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  scale_fill_manual(values = histology_colors) +
  theme_classic() +
  xlab("") + 
  ylab("L62 Residuals") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(fill = "none") +
  ggtitle("L62 (PCGA1)")
dev.off()

########## Module L70
pdf("./Figures/Fig6E_PCGA1_ModuleL70_Boxplot.pdf", height = 5, width = 4)
ggplot(pcga1.br.df, aes(x = Dysplasia_Grade, y = residuals_L70, fill = Dysplasia_Grade)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  scale_fill_manual(values = histology_colors) +
  theme_classic() +
  xlab("") + 
  ylab("L70 Expression") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(fill = "none") +
  ggtitle("L70 (PCGA1)")
dev.off()

##############################################################################################################
##################################### Figure 6F

############################## Modules
t <- read.table("../Figure5/072324_Biopsies+Brushes_noIonocytes_decontXcounts_zCellType_proteinCodingGenes_module_features_K40_L70.csv", sep = "\t", header = TRUE)
t2 <- as.list(t)
t3 <- lapply(t2, function(z){ z[!is.na(z) & z != ""]})
names(t3) <- gsub("V", "L", names(t3))
# basal_modules <- t3[which(names(t3) %in% c("L33", "L62", "L70"))]
basal_modules <- t3[which(names(t3) %in% c("L2", "L10", "L12", "L13", "L14", "L15", "L21", "L33", "L64", "L59", "L62", "L66", "L70"))]

############################## Read in data
histology_colors <- readRDS("../histology_colors.rds")
names(histology_colors)[which(names(histology_colors) == "Moderate Dysplasia")] <- "Dysplasia"
names(histology_colors)[2] <- "Normal"

resid_br <- readRDS("./Objects_Results/Fig6F_PCGA2_BronchialBrush_Residuals.rds")

### Run GSVA
# gsvaparam <- gsvaParam(exprData = resid_br, geneSets = basal_modules, kcdf = "Gaussian", absRanking = FALSE)
# br.gsva.scores <- gsva(param  = gsvaparam, verbose=TRUE)

### Save RDS files
# saveRDS(br.gsva.scores, "./Objects_Results/Fig6F_PCGA2_GSVA_Scores.rds")

### Read in GSVA objects
br.gsva.scores <- readRDS("./Objects_Results/Fig6F_PCGA2_GSVA_Scores.rds")

### Bulk annotations
bulk.se <- readRDS("./PCGA2_bulk_SE_object.rds")
rownames(bulk.se) <- bulk.se@rowRanges@elementMetadata@listData[["external_gene_name"]]
br.bulk <- bulk.se[,which(bulk.se$tissue_type == "Bronchial Brush")]
br.bulk.anndata <- as.data.frame(colData(br.bulk))

br.bulk <- cbind(br.bulk.anndata, t(br.gsva.scores))

### Edit metadata
br.bulk$Worst_Histology_Same_TimePoint_Clean <- gsub("Early squamous metaplasia", "Metaplasia", br.bulk$Worst_Histology_Same_TimePoint_Clean)
br.bulk$Worst_Histology_Same_TimePoint_Clean <- gsub("Basal Cell Hyperplasia", "Hyperplasia", br.bulk$Worst_Histology_Same_TimePoint_Clean)
br.bulk$Worst_Histology_Same_TimePoint_Clean <- gsub("Inflammation", "Normal", br.bulk$Worst_Histology_Same_TimePoint_Clean)
br.bulk$Worst_Histology_Same_TimePoint_Clean <- gsub("Mild Dysplasia", "Dysplasia", br.bulk$Worst_Histology_Same_TimePoint_Clean)
br.bulk$Worst_Histology_Same_TimePoint_Clean <- gsub("Moderate Dysplasia", "Dysplasia", br.bulk$Worst_Histology_Same_TimePoint_Clean)
br.bulk$Worst_Histology_Same_TimePoint_Clean <- gsub("Severe Dysplasia", "Dysplasia", br.bulk$Worst_Histology_Same_TimePoint_Clean)
br.bulk <- br.bulk[-which(br.bulk$Worst_Histology_Same_TimePoint_Clean == "Unknown"),]

br.bulk$Worst_Histology_Same_TimePoint_Clean <- factor(br.bulk$Worst_Histology_Same_TimePoint_Clean,
                                                       levels = c('Normal', 'Hyperplasia', 'Metaplasia', 'Dysplasia'))

br.bulk$Worst_Histology_Same_TimePoint_Clean2 <- as.numeric(br.bulk$Worst_Histology_Same_TimePoint_Clean)

############################## Modeling

lmRes_L33_br <- lme(L33 ~ Worst_Histology_Same_TimePoint_Clean2 + L70 + smoking.status.prediction, random = ~ 1|PCGA02.0_PatientID, data = br.bulk, na.action = na.omit)
summary(lmRes_L33_br)
saveRDS(lmRes_L33_br, "./Objects_Results/Fig6F_PCGA2_ModuleL33_linearModelResults.rds")

lmRes_L62_br <- lme(L62 ~ Worst_Histology_Same_TimePoint_Clean2 + L70 + smoking.status.prediction, random = ~ 1|PCGA02.0_PatientID, data = br.bulk, na.action = na.omit)
summary(lmRes_L62_br)
saveRDS(lmRes_L62_br, "./Objects_Results/Fig6F_PCGA2_ModuleL62_linearModelResults.rds")

lmRes_L70_br <- lme(L70 ~ Worst_Histology_Same_TimePoint_Clean2 + smoking.status.prediction, random = ~ 1|PCGA02.0_PatientID, data = br.bulk, na.action = na.omit)
summary(lmRes_L70_br)
saveRDS(lmRes_L70_br, "./Objects_Results/Fig6F_PCGA2_ModuleL70_linearModelResults.rds")

############# Save results from all modules
# Initialize a list to store results
modules <- c("L2", "L10", "L12", "L13", "L14", "L15", "L21", "L33", "L64", "L59", "L62", "L66")
model_results <- list()

for (mod in modules) {
  # Build the formula dynamically
  formula_str <- paste0(mod, " ~ Worst_Histology_Same_TimePoint_Clean2 + L70 + smoking.status.prediction")
  formula_obj <- as.formula(formula_str)
  
  # Fit the model
  model <- lme(formula_obj, random = ~ 1 | PCGA02.0_PatientID, data = br.bulk, na.action = na.omit)
  
  # Store the summary in a named list
  model_results[[mod]] <- summary(model)
}

model_results[["L70"]] <- summary(lmRes_L70_br)

saveRDS(model_results, "./Objects_Results/Fig6F_PCGA2_All6DModules_linearModelResults.rds")

############################## Plotting
model <- lme(L33 ~ L70 + smoking.status.prediction, random = ~ 1|PCGA02.0_PatientID, data = br.bulk, na.action = na.omit)
residuals_L33 <- residuals(model)
br.bulk$residuals_L33 <- residuals_L33

model <- lme(L62 ~ L70 + smoking.status.prediction, random = ~ 1|PCGA02.0_PatientID, data = br.bulk, na.action = na.omit)
residuals_L62 <- residuals(model)
br.bulk$residuals_L62 <- residuals_L62

model <- lme(L70 ~ smoking.status.prediction, random = ~ 1|PCGA02.0_PatientID, data = br.bulk, na.action = na.omit)
residuals_L70 <- residuals(model)
br.bulk$residuals_L70 <- residuals_L70

########## Module 33
pdf("./Figures/Fig6F_PCGA2_ModuleL33_Boxplot.pdf", height = 5, width = 4)
ggplot(br.bulk, aes(x = Worst_Histology_Same_TimePoint_Clean, y = residuals_L33, fill = Worst_Histology_Same_TimePoint_Clean)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  scale_fill_manual(values = histology_colors) +
  # scale_color_manual(values = c("red", "black"), labels = c("Br", "None")) +
  theme_classic() +
  xlab("") + 
  ylab("L33 Residuals") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(fill = "none", color = "none") +
  ggtitle("L33 (PCGA2)")
dev.off()

########## Module 62
pdf("./Figures/Fig6F_PCGA2_ModuleL62_Boxplot.pdf", height = 5, width = 4)
ggplot(br.bulk, aes(x = Worst_Histology_Same_TimePoint_Clean, y = residuals_L62, fill = Worst_Histology_Same_TimePoint_Clean)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  scale_fill_manual(values = histology_colors) +
  # scale_color_manual(values = c("red", "black"), labels = c("Br", "None")) +
  theme_classic() +
  xlab("") + 
  ylab("L62 Residuals") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(fill = "none", color = "none") +
  ggtitle("L62 (PCGA2)")
dev.off()

########## Module 70
pdf("./Figures/Fig6F_PCGA2_ModuleL70_Boxplot.pdf", height = 5, width = 4)
ggplot(br.bulk, aes(x = Worst_Histology_Same_TimePoint_Clean, y = residuals_L70, fill = Worst_Histology_Same_TimePoint_Clean)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  scale_fill_manual(values = histology_colors) +
  # scale_color_manual(values = c("red", "black"), labels = c("Br", "None")) +
  theme_classic() +
  xlab("") + 
  ylab("L70 Expression") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(fill = "none", color = "none") +
  ggtitle("L70 (PCGA2)")
dev.off()








##############################################################################################################
##################################### Figure 6D/E/F: all sig modules in PCGA1/2

PCGA1 <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/Fig4-6_Paper_Scripts_Results/Figure6/Objects_Results/Fig6E_PCGA1_All6DModules_linearModelResults.rds")
PCGA2 <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/Fig4-6_Paper_Scripts_Results/Figure6/Objects_Results/Fig6F_PCGA2_All6DModules_linearModelResults.rds")

res_table <- as.data.frame(matrix(0, nrow = 13, ncol = 5))

for (i in 1:13) {
  res_table[i,1] <- names(PCGA1)[i]
  res_table[i,2] <- PCGA1[[i]]$tTable[2,1]
  res_table[i,3] <- PCGA1[[i]]$tTable[2,5]
  res_table[i,4] <- PCGA2[[i]]$tTable[2,1]
  res_table[i,5] <- PCGA2[[i]]$tTable[2,5]
}

colnames(res_table) <- c("Module", "PCGA1_Coeff", "PCGA1_PValue", "PCGA2_Coeff", "PCGA2_PValue")





