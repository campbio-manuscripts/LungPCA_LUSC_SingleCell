library(celda)
library(singleCellTK)
library(SingleCellExperiment)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(GSVA)
library(RColorBrewer)
library(nlme)
library(ComplexHeatmap)
library(ggforce)
library(cowplot)
library(scales)
library(enrichR)
library(ggvenn)
library(ggtext)
library(sccomp)
library(openxlsx)
library(magrittr)
library(loo)
library(cmdstanr)
library(rlang)
library(parallel)
library(lifecycle)
library(readr)
library(glmmTMB)

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/")

sce_full <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/AllBiopsyCell_SCE_Final.rds")
sce <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/NonimmuneBiopsyCell_SCE_Final.rds")
celltype_colors <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/celltype_colors.rds") 
hist_cols <- readRDS("histology_colors.rds")
smoke_cols <- readRDS("smoking_colors.rds")
f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")), 
                                    colData = colData(altExp(sce)))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")

m <- match(colnames(sce), colnames(sce_full))
sce$CellTypeAllCells <- altExp(sce)$CellTypeAllCells <- sce_full$CellType[m]
sce$HistologyFull <- sce$Histology
sce$Histology <- as.character(sce$Histology)
sce$Histology[which(sce$Histology == "Denuded Bronchial Mucosa")] <- "Normal"
sce$Histology[which(sce$Histology == "Inflammation")] <- "Normal"
sce$Histology <- factor(sce$Histology, 
                        levels = c("Normal", "Basal Cell Hyperplasia", "Metaplasia", "Mild Dysplasia",
                                   "Moderate Dysplasia", "Severe Dysplasia", "CIS", "LUSC"),
                        ordered = TRUE)

sce$celda_cell_cluster_allcells <- sce_full$celda_cell_cluster[m]
sce$celda_cell_cluster <- celdaClusters(sce)
sce$CellType_Cluster <- sce$CellType
sce$CellType_Cluster <- gsub(" Cells","", sce$CellType_Cluster)
sce$CellType_Cluster <- gsub("Fibroblasts", "Fibroblast", sce$CellType_Cluster)
sce$CellType_Cluster <- paste(sce$CellType_Cluster, paste0("K",sce$celda_cell_cluster))
sce$CellType_Cluster <- factor(sce$CellType_Cluster,
                               levels = c(paste0("Basal K",1:9), "KRT5+/SCGB1A1+ K10", "KRT5+/MUC5B+ K11", "Perigoblet K12",
                                          paste0("Goblet K",13:17), paste0("Club K",18:22), "Mucous SMG K23", "Serous SMG K24",
                                          paste0("Ciliated K",25:30), paste0("Airway Smooth Muscle K",31), paste0("Fibroblast K",32:33), paste0("Endothelial K",34:35)),
                               ordered = TRUE)
sce$HistologyCondense <- as.character(sce$Histology)
sce$HistologyCondense[which(sce$Histology %in% c("Normal","Denuded Bronchial Mucosa","Inflammation"))] <- "NAD"
sce$HistologyCondense[which(sce$Histology %in% c("Basal Cell Hyperplasia","Metaplasia","Mild Dysplasia","Moderate Dysplasia","Severe Dysplasia"))] <- "PML"
sce$HistologyCondense[which(sce$Histology %in% c("CIS","LUSC"))] <- "CIS/LUSC"
sce$HistologyCondense <- factor(sce$HistologyCondense, levels = c("NAD","PML","CIS/LUSC"), ordered = TRUE)

## Figure 2A: Cell Type UMAP
sce$CellType_ClusterRange <- as.character(sce$CellType_Cluster)
k1.5 <- which(sce$celda_cell_cluster %in% c(1:5)); sce$CellType_ClusterRange[k1.5] <- "KRT5+ Basal Cells (K1-K5)"
k6 <- which(sce$celda_cell_cluster %in% c(6)); sce$CellType_ClusterRange[k6] <- "KRT5+ Basal Cells (K6)"
k7.9 <- which(sce$celda_cell_cluster %in% c(7:9)); sce$CellType_ClusterRange[k7.9] <- "KRT5+ Basal Cells (K7-K9)"
k10 <- which(sce$celda_cell_cluster %in% c(10)); sce$CellType_ClusterRange[k10] <- "KRT5+/SCGB1A1+ Cells (K10)"
k11 <- which(sce$celda_cell_cluster %in% c(11)); sce$CellType_ClusterRange[k11] <- "KRT5+/MUC5B+ Cells (K11)"
k12 <- which(sce$celda_cell_cluster %in% c(12)); sce$CellType_ClusterRange[k12] <- "CEACAM5+ Perigoblet Cells (K12)"
k13.17 <- which(sce$celda_cell_cluster %in% c(13:17)); sce$CellType_ClusterRange[k13.17] <- "MUC5AC+ Goblet Cells (K13-K17)"
k18.22 <- which(sce$celda_cell_cluster %in% c(18:22)); sce$CellType_ClusterRange[k18.22] <- "SCGB1A1+ Club Cells (K18-K22)"
k23 <- which(sce$celda_cell_cluster %in% c(23)); sce$CellType_ClusterRange[k23] <- "MUC5B+ Mucous SMG Cells (K23)"
k24 <- which(sce$celda_cell_cluster %in% c(24)); sce$CellType_ClusterRange[k24] <- "PRB1+ Serous SMG Cells (K24)"
k25.30 <- which(sce$celda_cell_cluster %in% c(25:30)); sce$CellType_ClusterRange[k25.30] <- "FOXJ1+ Ciliated Cells (K25-K30)"
k31 <- which(sce$celda_cell_cluster %in% c(31)); sce$CellType_ClusterRange[k31] <- "ACTA2+ Airway Smooth Muscle Cells (K31)"
k32.33 <- which(sce$celda_cell_cluster %in% c(32:33)); sce$CellType_ClusterRange[k32.33] <- "COL1A1+ Fibroblasts (K32-K33)"
k34.35 <- which(sce$celda_cell_cluster %in% c(34:35)); sce$CellType_ClusterRange[k34.35] <- "VWF+ Endothelial Cells (K34-K35)"
altExp(sce)$CellType_ClusterRange <- sce$CellType_ClusterRange

celltype_rangecols <- rep(NA, length(unique(sce$CellType_ClusterRange))); names(celltype_rangecols) <- unique(sce$CellType_ClusterRange)
celltype_rangecols["KRT5+ Basal Cells (K1-K5)"] <- celltype_colors[["Basal Cells"]]; celltype_rangecols["KRT5+ Basal Cells (K6)"] <- celltype_colors[["Basal Cells"]]; 
celltype_rangecols["KRT5+ Basal Cells (K7-K9)"] <- celltype_colors[["Basal Cells"]]; celltype_rangecols["KRT5+/SCGB1A1+ Cells (K10)"] <- celltype_colors[["KRT5+/SCGB1A1+ Cells"]];  
celltype_rangecols["KRT5+/MUC5B+ Cells (K11)"] <- celltype_colors[["KRT5+/MUC5B+ Cells"]]; celltype_rangecols["CEACAM5+ Perigoblet Cells (K12)"] <- celltype_colors[["Perigoblet Cells"]]; 
celltype_rangecols["MUC5AC+ Goblet Cells (K13-K17)"] <- celltype_colors[["Goblet Cells"]]; celltype_rangecols["SCGB1A1+ Club Cells (K18-K22)"] <- celltype_colors[["Club Cells"]]; 
celltype_rangecols["MUC5B+ Mucous SMG Cells (K23)"] <- celltype_colors[["Mucous SMG Cells"]]; celltype_rangecols["PRB1+ Serous SMG Cells (K24)"] <- celltype_colors[["Serous SMG Cells"]]; 
celltype_rangecols["FOXJ1+ Ciliated Cells (K25-K30)"] <- celltype_colors[["Ciliated Cells"]]; celltype_rangecols["ACTA2+ Airway Smooth Muscle Cells (K31)"] <- celltype_colors[["Airway Smooth Muscle Cells"]];  
celltype_rangecols["COL1A1+ Fibroblasts (K32-K33)"] <- celltype_colors[["Fibroblasts"]]; celltype_rangecols["VWF+ Endothelial Cells (K34-K35)"] <- celltype_colors[["Endothelial Cells"]]  

celltype_umap <- plotSCEDimReduceColData(altExp(sce), colorBy = "CellType_ClusterRange", reducedDimName = "celda_UMAP", dotSize = 0.25, xlab = "UMAP1", ylab = "UMAP2")
celltype_umap_altered <- celltype_umap + theme_classic() + scale_color_manual(values = celltype_rangecols) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), legend.position = "none")
pdf("Nonimmune_CellType_UMAP_wMarkers.pdf", width = 3, height = 3)
celltype_umap_altered
dev.off()
  
## Figure 2B: Marker Plots
markers <- c("KRT5", "KRT6A", "SCGB1A1", "MUC5B", "CEACAM5", "MUC5AC", "PRB1", "FOXJ1", "ACTA2", "COL1A1", "VWF")
modules <- read.table("070924_Nonimmune_celda_K35_L73_CeldaEdit.tsv", sep = "\t", header = TRUE)
modules <- lapply(as.list(modules), function(x) {
  x <- x[x != ""]
  x
})
modules_list <- as.list(modules)
markers_grep <- paste0(markers,"$")
marker_mods <- c()
for(i in markers) {
  marker_mods[i] <- grep(i, modules_list)
}
m <- match(paste0("L",marker_mods), rownames(sce_modular))
scaled_expr <- t(assay(sce_modular,"LogNormalize")[m,])
scaled_expr <- scale(scaled_expr)
colnames(scaled_expr) <- paste0("L",marker_mods," (",markers,")")
df <- cbind(data.frame(CellType = sce_modular$CellType), scaled_expr)
df_markersummary <- df %>% group_by(CellType) %>% summarize_all(mean)
df_markersummary[10,5] <- 4 ## Mucous SMG cells
df_markersummary[12,8] <- 4 ## Serous SMG cells
df_markersummary_pivot <- df_markersummary %>% pivot_longer(!CellType, names_to = "Module", values_to = "Expression")
df_markersummary_pivot$CellType <- factor(df_markersummary_pivot$CellType,
                                          levels = rev(c("Basal Cells", "KRT5+/SCGB1A1+ Cells", "KRT5+/MUC5B+ Cells", "Perigoblet Cells",
                                                         "Goblet Cells", "Club Cells", "Mucous SMG Cells", "Serous SMG Cells", "Ciliated Cells", 
                                                         "Airway Smooth Muscle Cells", "Fibroblasts", "Endothelial Cells")),
                                          ordered = TRUE)
df_markersummary_pivot$Module <- factor(df_markersummary_pivot$Module, levels = colnames(scaled_expr), ordered = TRUE)
midpoint <- mean(df_markersummary_pivot$Expression)
marker_scalehm <- ggplot(df_markersummary_pivot, aes(x = Module, y = CellType, fill = Expression)) + geom_tile() + 
  theme_minimal() + ylab("") + xlab("") +
  scale_fill_gradient2(name = "Average\nScaled\nNormalized\nExpression", midpoint = midpoint,
                       low = "blue", mid = "white", high = "red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, -0.5), "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.ticks = element_blank())
marker_scalehm

celltype_colorbar_df <- data.frame(CellType = levels(df_markersummary_pivot$CellType), X = 1)
celltype_colorbar_df$CellType <- factor(celltype_colorbar_df$CellType,
                                        levels = levels(df_markersummary_pivot$CellType), ordered = TRUE)
celltype_colorbar_ggbar <- ggplot(celltype_colorbar_df, aes(x = X, y = CellType, fill = CellType)) + geom_bar(position="fill", stat="identity") +
  theme_minimal() + scale_fill_manual(values = celltype_colors) + xlab("") + ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_blank(), legend.position = "none")
pdf("Marker_Modules_HM_Nonimmune_wCellType.pdf", width = 6, height = 4)
plot_grid(celltype_colorbar_ggbar, marker_scalehm, align = "h", nrow = 1, rel_widths = c(1,2))
dev.off()

## Compare relative proportions of Histology (High- vs. Low-Grade) and Smoking Status (Current vs. Smoker) in each cell type
sce$HistologyBinary <- ifelse(sce$Histology %in% c("Severe Dysplasia", "CIS", "LUSC"), "High Grade", "Low Grade")
sce$HistologyBinary <- factor(sce$HistologyBinary, levels = unique(sce$HistologyBinary), ordered = FALSE)
sce$Smoking_Status <- factor(sce$Smoking_Status, levels = unique(sce$Smoking_Status), ordered = FALSE)
sce$PCGA02_site <- factor(sce$PCGA02_site, levels = c("Roswell", "UCL"), ordered = FALSE)
sce$KCluster <- factor(sce$KCluster, levels = paste0("K",1:35), ordered = FALSE)
sce$Subject <- paste0("P", sce$Subject)
sce$Sample <- paste0("S", sce$Sample)
colData(sce) <- within(colData(sce), HistologyBinary <- relevel(HistologyBinary, ref = "Low Grade"))
colData(sce) <- within(colData(sce), Smoking_Status <- relevel(Smoking_Status, ref = "Former"))
colData(sce) <- within(colData(sce), PCGA02_site <- relevel(PCGA02_site, ref = "Roswell"))
colData(sce) <- within(colData(sce), KCluster <- relevel(KCluster, ref = "K1"))
sce_basal <- sce[, which(sce$KCluster %in% paste0("K",1:11))]
sce_basal$KCluster <- factor(sce_basal$KCluster, levels = paste0("K",1:11), ordered = FALSE)

cellshifts_binaryhistology_sccomp <- function(sce_object) {
  cat("Calculating Full Model.\n")
  model_with_factor_association <- 
    sce_object |>
    sccomp_estimate( 
      formula_composition = ~ 1 + HistologyBinary + Smoking_Status + PCGA02_site, 
      .sample = Sample, 
      .cell_group = KCluster, 
      inference_method = "hmc",
      enable_loo = TRUE
    ) |> 
    sccomp_remove_outliers(cores = 4) |> 
    sccomp_test()
  plot_sum_histology_smoking <- model_with_factor_association |> 
    sccomp_boxplot(factor = "HistologyBinary")
  list(fullmodel = model_with_factor_association, fullmodel_plot = plot_sum_histology_smoking)
}
nonimmune_binaryhistology_cellshifts <- cellshifts_binaryhistology_sccomp(sce_object = sce)
nonimmune_model_res <- nonimmune_binaryhistology_cellshifts$fullmodel %>% arrange(parameter, KCluster) %>% select(KCluster:c_FDR)
nonimmune_model_res[is.na(nonimmune_model_res)] <- ""

## Classify samples based on basal cell repertoire
sce$BasalCluster <- sce$HGBasalPercent <- sce$HGBasalPercentMixed <- sce$BasalGroup <- sce$BasalGroupMixed <- sce$TotalBasalNumbers <- sce$TotalBasalNumbersMixed <- c()
sce$BasalCluster <- ""
sce$BasalCluster[which(sce$celda_cell_cluster %in% c(1:6))] <- "High Grade Basal Cells"
sce$BasalCluster[which(sce$celda_cell_cluster %in% c(7:9))] <- "Low Grade Basal Cells"
sce$BasalCluster <- factor(sce$BasalCluster, levels = c("Low Grade Basal Cells", "High Grade Basal Cells"), ordered = TRUE)
sce$BasalClusterMixed <- ""
sce$BasalClusterMixed[which(sce$celda_cell_cluster %in% c(1:6))] <- "High Grade Basal Cells"
sce$BasalClusterMixed[which(sce$celda_cell_cluster %in% c(7:11))] <- "Low Grade Basal Cells"
sce$BasalClusterMixed <- factor(sce$BasalClusterMixed, levels = c("Low Grade Basal Cells", "High Grade Basal Cells"), ordered = TRUE)

cd <- data.frame(colData(sce))
cd$CellType[which(cd$celda_cell_cluster %in% c(1:6))] <- "High Grade Basal Cells"
basalnumber_sample <- cd %>% group_by(Sample) %>% dplyr::count(CellType, .drop = FALSE)
basalnumber_sample <- merge(basalnumber_sample, unique(data.frame(Sample = cd$Sample, Histology = cd$Histology, SmokingStatus = cd$Smoking_Status, row.names = colnames(sce))), by = "Sample")
basalnumber_sample$LGBasalNumbers <- rep(0, times = nrow(basalnumber_sample))
basalnumber_sample$LGBasalNumbersMixed <- rep(0, times = nrow(basalnumber_sample))
basalnumber_sample$HGBasalNumbers <- rep(0, times = nrow(basalnumber_sample))
basalnumber_sample$TotalBasalNumbers <- rep(0, times = nrow(basalnumber_sample))
basalnumber_sample$TotalBasalNumbersMixed <- rep(0, times = nrow(basalnumber_sample))
for(i in 1:nrow(basalnumber_sample)) {
  high_grade_basal_num <- basalnumber_sample[which(basalnumber_sample$Sample == basalnumber_sample$Sample[i] & basalnumber_sample$CellType == "High Grade Basal Cells"),"n"]
  low_grade_basal_num <- basalnumber_sample[which(basalnumber_sample$Sample == basalnumber_sample$Sample[i] & basalnumber_sample$CellType == "Basal Cells"),"n"]
  low_grade_basal_num_mixed <- sum(basalnumber_sample[which(basalnumber_sample$Sample == basalnumber_sample$Sample[i] & basalnumber_sample$CellType %in% c("Basal Cells", "KRT5+/SCGB1A1+ Cells", "KRT5+/MUC5B+ Cells")),"n"])
  if(length(high_grade_basal_num) == 0) {
    high_grade_basal_num <- 0
  }
  if(length(low_grade_basal_num) == 0) {
    low_grade_basal_num <- 0
  }
  if(length(low_grade_basal_num_mixed) == 0) {
    low_grade_basal_num_mixed <- 0
  }
  basalnumber_sample$LGBasalNumbers[i] <- low_grade_basal_num
  basalnumber_sample$LGBasalNumbersMixed[i] <- low_grade_basal_num_mixed
  basalnumber_sample$HGBasalNumbers[i] <- high_grade_basal_num
  basalnumber_sample$TotalBasalNumbers[i] <- sum(low_grade_basal_num + high_grade_basal_num)
  basalnumber_sample$TotalBasalNumbersMixed[i] <- sum(low_grade_basal_num_mixed + high_grade_basal_num)
}
basalnumber_sample$HGBasalPercent <- basalnumber_sample$HGBasalNumbers/basalnumber_sample$TotalBasalNumbers
basalnumber_sample$HGBasalPercent[which(is.nan(basalnumber_sample$HGBasalPercent))] <- 0
basalnumber_sample$HGBasalPercentMixed <- basalnumber_sample$HGBasalNumbers/basalnumber_sample$TotalBasalNumbersMixed
basalnumber_sample$HGBasalPercentMixed[which(is.nan(basalnumber_sample$HGBasalPercentMixed))] <- 0

# Transfer Basal Group labels to Nonimmune SCE
sce$HGBasalPercent <- 0
sce$HGBasalPercentMixed <- 0
sce$TotalBasalNumbers <- 0
sce$TotalBasalNumbersMixed <- 0
for(i in unique(basalnumber_sample$Sample)) {
  samp_index <- which(basalnumber_sample$Sample == i)[1]
  sce_index <- which(sce$Sample == i)
  sce$HGBasalPercent[sce_index] <- basalnumber_sample[samp_index, "HGBasalPercent"]
  sce$HGBasalPercentMixed[sce_index] <- basalnumber_sample[samp_index, "HGBasalPercentMixed"]
  sce$TotalBasalNumbers[sce_index] <- basalnumber_sample[samp_index, "TotalBasalNumbers"]
  sce$TotalBasalNumbersMixed[sce_index] <- basalnumber_sample[samp_index, "TotalBasalNumbersMixed"]
}

sce$BasalGroup <- ""
sce$BasalGroup[which(sce$HGBasalPercent < 0.2)] <- "Low Grade Basal Sample"
sce$BasalGroup[which(sce$HGBasalPercent > 0.2 & sce$HGBasalPercent < 0.8)] <- "Mixed Grade Basal Sample"
sce$BasalGroup[which(sce$HGBasalPercent > 0.8)] <- "High Grade Basal Sample"
sce$BasalGroup[which(sce$TotalBasalNumbers < 10)] <- "Too Few Basal Cells"
sce$BasalGroupMixed <- ""
sce$BasalGroupMixed[which(sce$HGBasalPercentMixed < 0.2)] <- "Low Grade Basal Sample"
sce$BasalGroupMixed[which(sce$HGBasalPercentMixed > 0.2 & sce$HGBasalPercentMixed < 0.8)] <- "Mixed Grade Basal Sample"
sce$BasalGroupMixed[which(sce$HGBasalPercentMixed > 0.8)] <- "High Grade Basal Sample"
sce$BasalGroupMixed[which(sce$TotalBasalNumbersMixed < 10)] <- "Too Few Basal Cells"

# Transfer Basal Group to All Cell SCE
sce_full$HGBasalPercent <- sce_full$HGBasalPercentMixed <- sce_full$TotalBasalNumbers <- sce_full$TotalBasalNumbersMixed <- 0
sce_full$BasalGroup <- sce_full$BasalGroupMixed <- ""
for(i in unique(sce_full$Sample)) {
  if(!i %in% unique(sce$Sample)) {
    noepi <- which(sce_full$Sample == i)
    sce_full$HGBasalPercent[noepi] <- sce_full$HGBasalPercentMixed[noepi] <- sce_full$TotalBasalNumbers[noepi] <- sce_full$TotalBasalNumbersMixed[noepi] <- 0
    sce_full$BasalGroup[noepi] <- sce_full$BasalGroupMixed[noepi] <- "No Epithelial Cells"
  } else {
    s_total <- which(sce_full$Sample == i)
    s_epi <- which(sce$Sample == i)
    sce_full$HGBasalPercent[s_total] <- unique(sce$HGBasalPercent[s_epi])
    sce_full$HGBasalPercentMixed[s_total] <- unique(sce$HGBasalPercentMixed[s_epi])
    sce_full$TotalBasalNumbers[s_total] <- unique(sce$TotalBasalNumbers[s_epi])
    sce_full$TotalBasalNumbersMixed[s_total] <- unique(sce$TotalBasalNumbersMixed[s_epi])
    sce_full$BasalGroup[s_total] <- unique(sce$BasalGroup[s_epi])
    sce_full$BasalGroupMixed[s_total] <- unique(sce$BasalGroupMixed[s_epi])
  }
}

## Figure 2C: Connect Basal Cell Cells to Histology
cluster_celltype_basalgroup_colorbar_df <- unique(data.frame(CellType = sce$CellType, Cluster = sce$KCluster, BasalCluster = sce$BasalClusterMixed, X = 1))
cluster_celltype_basalgroup_colorbar_df %<>% subset(!is.na(BasalCluster))
cluster_celltype_basalgroup_colorbar_df$Cluster <- factor(cluster_celltype_basalgroup_colorbar_df$Cluster,
                                                          levels = paste0("K",11:1), ordered = TRUE)
cluster_celltype_ggbar <- ggplot(cluster_celltype_basalgroup_colorbar_df, aes(x = X, y = Cluster, fill = CellType)) + geom_bar(position="fill", stat="identity") +
  theme_minimal() + scale_fill_manual(name = "Cell Type", values = celltype_colors) + xlab("Cell Type") + ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom", axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  guides(fill = guide_legend(nrow = 3, reverse = TRUE))

cluster_basalcluster_ggbar <- ggplot(cluster_celltype_basalgroup_colorbar_df, aes(x = X, y = Cluster, fill = BasalCluster)) + geom_bar(position="fill", stat="identity") +
  theme_minimal() + scale_fill_manual(name = "Basal Cell Grade", values = celltype_colors) + xlab("Basal Cell Grade") + ylab("Basal Cell Cluster") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom", axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  guides(fill = guide_legend(nrow = 2, reverse = TRUE))

df_hist <- data.frame(table(sce$KCluster, sce$Histology))
df_hist$Var1 <- factor(df_hist$Var1, levels = paste0("K",35:1), ordered = TRUE)
df_hist$Var2 <- factor(df_hist$Var2, levels = rev(levels(sce$Histology)), ordered = TRUE)
df_hist_basalonly <- df_hist[which(df_hist$Var1 %in% paste0("K",1:11)),]
histology_sbp_basalonly <- ggplot(df_hist_basalonly, aes(fill=Var2, x=Freq, y=Var1)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() + 
  ylab("") + xlab("% Histology") + ggtitle("") +
  scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), expand = expansion(add = c(0,0.02))) + scale_fill_manual(name = "Histology", values = hist_cols) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
  guides(fill = guide_legend(nrow = 5, reverse = TRUE))

df_smoke <- data.frame(table(sce$KCluster, sce$Smoking_Status))
df_smoke$Var1 <- factor(df_smoke$Var1, levels = paste0("K",35:1), ordered = TRUE)
df_smoke$Var2 <- factor(df_smoke$Var2, levels = rev(levels(sce$Smoking_Status)), ordered = TRUE)
df_smoke_basalonly <- df_smoke[which(df_smoke$Var1 %in% paste0("K",1:11)),]
smoke_sbp_basalonly <- ggplot(df_smoke_basalonly, aes(fill=Var2, x=Freq, y=Var1)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() + 
  ylab("") + xlab("% Smoking") + ggtitle("") +
  scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), expand = expansion(add = c(0,0.02))) + scale_fill_manual(name = "Smoking Status", values = smoke_cols) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
  guides(fill = guide_legend(nrow = 2, reverse = TRUE))

stats_values <- nonimmune_model_res %>% filter(KCluster %in% paste0("K",1:11)) %>% filter(factor %in% c("HistologyBinary", "Smoking_Status")) # If doing only Basal Cells
stats_values$c_FDR[stats_values$c_FDR == 0] <- 2e-16
stats_values$c_effect[which(stats_values$c_effect < -3)] <- -3 
stats_values$c_effect[which(stats_values$c_effect > 3)] <- 3 
stats_values$logp <- -log10(stats_values$c_FDR)
stats_values$sig <- NA
for (i in 1:length(stats_values$c_FDR)) {
  if (is.na(stats_values$c_FDR[i])) 
    stats_values$sig[i] <- "" 
  else if (stats_values$c_FDR[i] > 0.05)
    stats_values$sig[i] <- ""
  else if (stats_values$c_FDR[i] > 0.01) 
    stats_values$sig[i] <- "*" 
  else if (stats_values$c_FDR[i] > 0.001) 
    stats_values$sig[i] <- "**"  
  else stats_values$sig[i] <- "***"
}

stats_values$parameter <- gsub("HistologyBinaryHigh Grade", "Severe Dysplasia - LUSC vs. Normal - Moderate Dysplasia", stats_values$parameter)
stats_values$parameter <- gsub("Smoking_StatusCurrent", "Current vs. Former Smoker", stats_values$parameter)
stats_values$parameter <- factor(stats_values$parameter,
                                 levels = c("Severe Dysplasia - LUSC vs. Normal - Moderate Dysplasia", "Current vs. Former Smoker"),
                                 ordered = TRUE)
stats_values$KCluster <- factor(stats_values$KCluster,
                                levels = paste0("K",35:1), # Change to 11 if just plotting basal cells as below.
                                ordered = TRUE)

colors <- colorRampPalette(c("#053061", "#053061", "#053061", "#2166AC", "#4393C3", 
                             "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#67001F", 
                             "#67001F", "#67001F"))(200)
cellshift_cluster_histology_bubble <- ggplot(stats_values, aes(x = parameter, y = KCluster, label = sig, color = c_effect, 
                                             size = ifelse(is.na(logp), 2, logp), shape = ifelse(is.na(c_effect), "Missing", "Present"))) + 
  geom_point(pch=1, colour = "black", stroke = 2) +
  geom_point(stat = "identity", stroke = 2, na.rm = T, shape = 19) +  
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5),
        axis.ticks = element_blank()) +
  labs(x = "", y = "", color = "Effect Size", size = "-log10(FDR)") + 
  guides(size = "legend", color = "colorbar") +
  geom_text(color = "white", size = 3, vjust = 0.9) +
  scale_x_discrete(labels = c("Severe Dysplasia - LUSC\nvs.\nNormal - Moderate Dysplasia", "Current vs.\nFormer Smoker")) +
  scale_shape_manual(values = c(Missing = 4, Present = 19)) +
  scale_color_gradientn(colours = colors, limits = c(-3, 3), breaks = c(-3, 0, 3)) 

# Combine boxplots and sccomp figure
pdf("Basal_Cell_Cluster_Hist_Barplots+SCComp_wSite.pdf", width = 20, height = 6)
plot_grid(plotlist = list(cluster_basalcluster_ggbar, cluster_celltype_ggbar, histology_sbp_basalonly, smoke_sbp_basalonly, cellshift_cluster_histology_bubble), align = "hv", nrow = 1, rel_widths = c(1,1,2,2,4))
dev.off()

## Figure 2D: Basal Cluster UMAP
umap_coords <- reducedDim(altExp(sce), "celda_UMAP")
colnames(umap_coords) <- c("UMAP1", "UMAP2")
basalcluster_df <- cbind(data.frame(BasalCluster = sce$BasalCluster), umap_coords)
basalcluster_umap <- ggplot(basalcluster_df, aes(x = UMAP1, y = UMAP2, color = BasalCluster)) + geom_point(size = 0.25) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",
                          axis.ticks = element_blank()) + scale_color_manual(name = "Basal Cluster", breaks = c("Low Grade Basal Cells", "High Grade Basal Cells"), values = basalcluster_cols, na.value = "gray80") 
pdf("Basal_Cluster_UMAP_wMixed.pdf", width = 3, height = 3)
basalcluster_umap 
dev.off()

## LGB vs HGB modules - classifying double positive "Mixed" clusters 10-11 as basal cells
cellprobs <- data.frame(t(assay(sce_modular,"LogNormalize")))
cellprobs$BasalCluster <- sce$BasalClusterMixed
cellprobs$Sample <- sce$Sample
cellprobs$PCGA02_site <- sce$PCGA02_site
cellprobs$Smoking_Status <- sce$Smoking_Status
mean_cellprobs <- cellprobs %>% select(-c(Smoking_Status, Sample, PCGA02_site)) %>% group_by(BasalCluster) %>% summarize_all(mean)
mean_cellprobs <- data.frame(t(mean_cellprobs))
colnames(mean_cellprobs) <- mean_cellprobs["BasalCluster",]
mean_cellprobs <- mean_cellprobs[2:74,1:2] # remove non-basal cells
mean_cellprobs$`Low Grade Basal Cells` <- as.numeric(mean_cellprobs$`Low Grade Basal Cells`); mean_cellprobs$`High Grade Basal Cells` <- as.numeric(mean_cellprobs$`High Grade Basal Cells`)
mean_cellprobs$Log2FC <- log2(mean_cellprobs$`High Grade Basal Cells`/mean_cellprobs$`Low Grade Basal Cells`)
mean_cellprobs$Module <- rownames(mean_cellprobs)
colnames(mean_cellprobs) <- c("Mean_LGB", "Mean_HGB", "Log2FC", "Module")

# Calculate P-values for module differences
cellprobs_basalonly <- cellprobs
cellprobs_basalonly$BasalCluster <- factor(cellprobs_basalonly$BasalCluster, levels = unique(cellprobs_basalonly$BasalCluster), ordered = FALSE)
cellprobs_basalonly$Smoking_Status <- factor(cellprobs_basalonly$Smoking_Status, levels = unique(cellprobs_basalonly$Smoking_Status), ordered = FALSE)
cellprobs_basalonly$PCGA02_site <- factor(cellprobs_basalonly$PCGA02_site, levels = unique(cellprobs_basalonly$PCGA02_site), ordered = FALSE)
cellprobs_basalonly <- within(cellprobs_basalonly, BasalCluster <- relevel(BasalCluster, ref = "Low Grade Basal Cells"))
cellprobs_basalonly <- within(cellprobs_basalonly, Smoking_Status <- relevel(Smoking_Status, ref = "Former"))
cellprobs_basalonly <- within(cellprobs_basalonly, PCGA02_site <- relevel(PCGA02_site, ref = "Roswell"))
cellprobs_basalonly <- cellprobs_basalonly[which(cellprobs_basalonly$BasalCluster %in% c("Low Grade Basal Cells", "High Grade Basal Cells")),]
model_p <- data.frame(Module = names(cellprobs_basalonly)[1:73], row.names = names(cellprobs_basalonly)[1:73],
                      CellType_Beta = rep(0, 73), CellType_P = rep(0, 73), Smoking_Beta = rep(0, 73), Smoking_P = rep(0, 73))
for(i in rownames(model_p)) {
  model <- lme(as.formula(paste0(i," ~ BasalCluster + Smoking_Status + PCGA02_site")),
               random = ~ 1|Sample,
               data = cellprobs_basalonly,
               na.action = na.omit)
  tTab <- summary(model)$tTable
  model_p[i,"CellType_Beta"] <- tTab["BasalClusterHigh Grade Basal Cells","Value"]
  model_p[i,"CellType_P"] <- tTab["BasalClusterHigh Grade Basal Cells","p-value"]
  model_p[i,"Smoking_Beta"] <- tTab["Smoking_StatusCurrent","Value"]
  model_p[i,"Smoking_P"] <- tTab["Smoking_StatusCurrent","p-value"]
}
model_p$CellType_Q <- p.adjust(model_p$CellType_P, method = "fdr")
model_p$Smoking_Q <- p.adjust(model_p$Smoking_P, method = "fdr")
rm(model, tTab)
volcano_df <- merge(mean_cellprobs, model_p)

## Figure 2E: Heatmap of LGB vs. HGB modules
f_basal <- cellprobs[which(celdaClusters(sce) %in% c(1:11)),]
significant_modules <- volcano_df$Module[which(volcano_df$Label == "Label")]
cell_probs_scale_mat <- t(scale(f_basal[,significant_modules]))
cell_probs_scale_mat2 <- cell_probs_scale_mat
cell_probs_scale_mat2[cell_probs_scale_mat2 > 2] <- 2; cell_probs_scale_mat2[cell_probs_scale_mat2 < -2] <- -2
colorScheme <- circlize::colorRamp2(c(min(cell_probs_scale_mat2),(max(cell_probs_scale_mat2) + min(cell_probs_scale_mat2))/2, max(cell_probs_scale_mat2)),
                                    c("blue", "white","red"))
sce_basal <- sce[, which(celdaClusters(sce) %in% c(1:11))]
cellAnnot <- data.frame("Basal Cell Grade" = sce_basal$BasalClusterMixed, "Cluster" = sce_basal$KCluster, "Sample" = sce_basal$Sample, "Smoking Status" = sce_basal$Smoking_Status, "Histology" = sce_basal$Histology, row.names = colnames(sce_basal), check.names = FALSE)
cellAnnot_noclust <- cellAnnot[,c(1,3:5)]
sample_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(sce_basal$Sample))); names(sample_cols) <- unique(cellAnnot$Sample)
basalcluster_cols <- celltype_colors[c("Low Grade Basal Cells", "High Grade Basal Cells")]
CellAnnotColor <- list("Basal Cell Grade" = basalcluster_cols, "Sample" = sample_cols, "Smoking Status" = smoke_cols, Histology = hist_cols) # see Figure 1.R for cluster_colors2, also it's saved
ca <- ComplexHeatmap::HeatmapAnnotation(df = cellAnnot_noclust, 
                                        col = CellAnnotColor, show_legend = c(TRUE, FALSE, TRUE, TRUE))
cell_probs_hm <- ComplexHeatmap::Heatmap(matrix = cell_probs_scale_mat2, col = colorScheme, row_title = "Module", column_title = "Cell",
                                         cluster_rows = TRUE, cluster_columns = TRUE, column_split = sce_basal$KCluster, 
                                         show_column_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, top_annotation = ca,
                                         heatmap_legend_param = list(title = "Cell\nProbability"))
ht <- draw(cell_probs_hm); col_order_list <- column_order(ht); row_order_list <- row_order(ht)
col_order <- c(col_order_list$K7, col_order_list$K8, col_order_list$K9, col_order_list$K10, col_order_list$K11, # LGB
               col_order_list$K1, col_order_list$K2, col_order_list$K3, col_order_list$K4, col_order_list$K5, col_order_list$K6) # HGB
cell_probs_scale_mat2 <- cell_probs_scale_mat2[rev(row_order_list),col_order]
cellAnnot <- cellAnnot[col_order,]
cellAnnot_noclust <- cellAnnot_noclust[col_order,]
ca <- ComplexHeatmap::HeatmapAnnotation(df = cellAnnot_noclust, 
                                         col = CellAnnotColor, show_legend = c(TRUE, FALSE, TRUE, TRUE))
cell_probs_hm <- ComplexHeatmap::Heatmap(matrix = cell_probs_scale_mat2, col = colorScheme, row_title = "", column_title = "",
                                         cluster_rows = FALSE, cluster_columns = FALSE, column_split = cellAnnot$Cluster, 
                                         show_column_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, top_annotation = ca,
                                         heatmap_legend_param = list(title = "Scaled\nNormalized\nExpression"))
pdf("LGB_HGB_Module_Heatmap_wMixed_wSite.pdf", width = 11, height = 7)
cell_probs_hm
dev.off()

## Figure 2F: Violins for Select LGB vs. HGB Modules
basal_kcluster_cols <- c(rep("red",6),rep("blue",5))
names(basal_kcluster_cols) <- paste0("K",1:11)
significant_modules_interesting <- c("L17","L23","L1","L32")
cellprobs_significant_interesting <- lapply(X = significant_modules_interesting, FUN = function(mod, BasalCluster = "KCluster") {
  cellprobs_sig <- cellprobs_basalonly[,c(mod, BasalCluster)]
  cellprobs_sig[[BasalCluster]] <- factor(cellprobs_sig[[BasalCluster]],
                                          levels = c(paste0("K",1:11)),
                                          ordered = TRUE)
  sigmod_violin <- ggplot(cellprobs_sig, aes_string(x = BasalCluster, y = mod, fill = BasalCluster)) + geom_violin(scale = "width") +
    theme_classic() + labs(title = mod, y = "Normalized Expression", x = "Basal Cell Cluster") + scale_fill_manual(values = basal_kcluster_cols) + ylim(c(0,8)) + 
    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_markdown(hjust = 0.5), legend.position = "none") 
  sigmod_violin
})
rename_plot <- function(plt, title) {
  plt <- plt + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
}
cellprobs_significant_interesting[[1]] <- rename_plot(plt = cellprobs_significant_interesting[[1]],
                                                      title = "L17 - CDKN2A, MAGE")
cellprobs_significant_interesting[[2]] <- rename_plot(plt = cellprobs_significant_interesting[[2]],
                                                      title = "L23 - MKI67, TOP2A")
cellprobs_significant_interesting[[3]] <- rename_plot(plt = cellprobs_significant_interesting[[3]],
                                                      title = "L1 - MHCII")
cellprobs_significant_interesting[[4]] <- rename_plot(plt = cellprobs_significant_interesting[[4]],
                                                      title = "L32 - BCAM, KRT15")
pdf("HGB_vs_LGB_Module_Violins_Interesting_wMixed.pdf", width = 5, height = 4)
plot_grid(plotlist = cellprobs_significant_interesting, nrow = 2, align = "hv") 
dev.off()

## Score Outside Modules in Single Cell Data
# Load modules
gs_direction <- readRDS("/restricted/projectnb/pcga/ANALYSIS_NEXTFLOW/cor_modules/combined_gene_set_symbols_wdirME.rds")
gs_direction$bx_purple_neg <- c()  # This directional module has 2 genes
merrick_progreg <- openxlsx::read.xlsx("MerrickD_Histology-PersistVsRegress_MeanHisto_AllGenes.xlsx", sheet = "PBD_vs_RBD_forR")
merrick_progreg <- merrick_progreg[which(merrick_progreg$FDR != "?"),]
merrick_progreg_progup <- merrick_progreg[which(merrick_progreg$Relative.Change > 0),]
merrick_progreg_progdown <- merrick_progreg[which(merrick_progreg$Relative.Change < 0),]
merrick_hist <- openxlsx::read.xlsx("MerrickD_Histology-PersistVsRegress_MeanHisto_AllGenes.xlsx", sheet = "HistCorr_forR")
merrick_hist_hgup <- merrick_hist[which(merrick_hist$`HS.Corr.r-value` > 0),]
merrick_hist_hgdown <- merrick_hist[which(merrick_hist$`HS.Corr.r-value` < 0),]
merrick_list <- list("Merrick_ProgReg_Up" = merrick_progreg_progup$Gene.Symbol, "Merrick_ProgReg_Down" = merrick_progreg_progdown$Gene.Symbol,
                     "Merrick_HG_Up" = merrick_hist_hgup$Gene.Symbol, "Merrick_HG_Down" = merrick_hist_hgdown$Gene.Symbol)
mascaux <- read.table("final_modules_wunique_Gene_Cluster_pairs.txt", header = TRUE)
colnames(mascaux) <- c("Probe.ID", "Gene.Name", "Module")
mascaux_list <- split(mascaux, ~ Module)
mascaux_list <- lapply(mascaux_list, function(x) {
  x <- x[["Gene.Name"]]
})
names(mascaux_list) <- paste0("Mascaux_M", names(mascaux_list))
mascaux_list$Mascaux_M0 <- c()
# Identify Mascaux Modules
mascaux_enrichr <- list()
mascaux_enrichr <- lapply(mascaux_list, function(x) {
  df <- enrichr(genes = x, databases = "MSigDB_Hallmark_2020")
  df$MSigDB_Hallmark_2020
})

janes_epistromagenes_pvr <- openxlsx::read.xlsx("Janes_EpiStromaGenes_PvR.xlsx", sheet = "PvR_GXN")
janes_epistroma_progup <- janes_epistromagenes_pvr$gene[which(janes_epistromagenes_pvr$fc > 0)]
janes_epistroma_progdown <- janes_epistromagenes_pvr$gene[which(janes_epistromagenes_pvr$fc < 0)]
janes_cin <- c("ACTL6A", "ELAVL1", "MAD2L1", "NEK2", "OIP5")
janes_immune_genes_progreg <- c("B2M","HLA-A","HLA-B","HLA-C","HLA-DQA1","HLA-DPA1","JAK2","IRF2","IRF1","IRF3","PVR","BCL6","JAK1","CD80","CD86","HHLA2","IRF4","TNFSF18","TNFSF4","IRF5","IRF7","CD40","CD70","IRF8","TNFSF14","TNFSF9","RAET1E","ULBP1","IKBKG","CD274","PDCD1LG2","KDR","IDO1","IKBKB","IL10RB","LGALS3","CD276","VTCN1","CHUK","NOS2","CANX","HSPA5","PDIA3","CALR")
janes_immune_genes_progup <- c("B2M","HLA-B","PVR","BCL6","JAK1","TNFSF18","CD40","CD70","IRF8","RAET1E","ULBP1","IKBKG","PDCD1LG2","KDR","IKBKB","LGALS3","CD276","VTCN1","CHUK","CANX","HSPA5","CALR")
janes_immune_genes_progdown <- setdiff(janes_immune_genes_progreg, janes_immune_genes_progup)
janes_list <- list("Janes_epistroma_progup" = janes_epistroma_progup, "Janes_epistroma_progdown" = janes_epistroma_progdown, "Janes_cin" = janes_cin,
                   "Janes_immune_genes_progup" = janes_immune_genes_progup, "Janes_immune_genes_progdown" = janes_immune_genes_progdown)
full_module_list <- c(gs_direction, merrick_list, mascaux_list, janes_list)
full_module_list <- lapply(full_module_list, unique)
full_module_list <- lapply(full_module_list, function(x) {
  x <- x[!is.na(x)]
})

sce_scoring <- SingleCellExperiment(assays = SimpleList("decontXcounts" = assay(sce_full, "decontXcounts"), "LogNormalize" = assay(sce_full, "LogNormalize")),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce_full), "celda_UMAP")),
                                    colData = colData(sce_full),
                                    rowData = rowData(sce_full))
# Implement expression filter
sce_scoring <- selectFeatures(sce_scoring, useAssay = "decontXcounts") 
full_module_list_filtered <- lapply(full_module_list, function(x) {
  x <- x[x %in% rowData(altExp(sce_scoring))$Gene_name]
})

altExp(sce_scoring) <- importGeneSetsFromList(altExp(sce_scoring), geneSetList = full_module_list_filtered, collectionName = "All_Mods", by = "Gene_name")
altExp(sce_scoring) <- runVAM(inSCE = altExp(sce_scoring),
                      geneSetCollectionName = "All_Mods",
                      useAssay = "LogNormalize")
vam_scores <- data.frame(reducedDim(altExp(sce), "VAM_All_Mods_CDF"))
saveRDS(vam_scores, "Outside_Modules_VAM_Scores_All_Biopsy_Cells.rds")

# Because VAM scores are 0-1 restricted, modeling them as a continuous variable (i.e., as a function of categorical variables) violates the assumptions of linear modeling. Thus, using a logistic regression model
vam_model_logistic <- function(scores) {
  scores_nonimmune <- scores[colnames(sce),]
  scores_nonimmune$CellType <- sce$CellType
  scores_nonimmune$Smoking_Status <- sce$Smoking_Status
  scores_nonimmune$Sample <- sce$Sample
  scores_nonimmune$KCluster <- sce$KCluster
  scores_nonimmune$PCGA02_site <- sce$PCGA02_site
  scores_nonimmune$CellType[which(scores_nonimmune$KCluster %in% paste0("K",1:6))] <- "High Grade Basal Cells"
  scores_nonimmune$CellType[which(scores_nonimmune$KCluster %in% paste0("K",7:11))] <- "Low Grade Basal Cells"
  scores_basal <- scores_nonimmune[scores_nonimmune$CellType %in% c("High Grade Basal Cells", "Low Grade Basal Cells"),]
  scores_basal$CellType <- factor(scores_basal$CellType, levels = c("Low Grade Basal Cells", "High Grade Basal Cells"), ordered = FALSE)
  scores_basal$Smoking_Status <- factor(scores_basal$Smoking_Status, levels = levels(sce$Smoking_Status), ordered = FALSE)
  scores_basal$PCGA02_site <- factor(scores_basal$PCGA02_site, levels = unique(sce$PCGA02_site), ordered = FALSE)
  scores_basal <- within(scores_basal, CellType <- relevel(CellType, ref = "Low Grade Basal Cells"))
  scores_basal <- within(scores_basal, Smoking_Status <- relevel(Smoking_Status, ref = "Former"))
  scores_basal <- within(scores_basal, PCGA02_site <- relevel(PCGA02_site, ref = "Roswell"))
  model_p <- data.frame(Module = names(full_module_list_filtered),
                        LGB_Mean = rep(0, 28),
                        HGB_Mean = rep(0, 28),
                        Score_Logit = rep(0, 28),
                        Score_P = rep(0, 28),
                        Smoke_Logit = rep(0, 28),
                        Smoke_P = rep(0, 28),
                        row.names = names(full_module_list_filtered))
  for(i in model_p$Module) {
    model_p[i, "LGB_Mean"] <- mean(scores_basal[which(scores_basal$CellType == "Low Grade Basal Cells"),i])
    model_p[i, "HGB_Mean"] <- mean(scores_basal[which(scores_basal$CellType == "High Grade Basal Cells"),i])
    model_terms <- paste0("CellType ~ ", i, " + Smoking_Status + PCGA02_site + (1|Sample)")
    model_terms <- as.formula(model_terms)
    lmres <- glmmTMB(model_terms,
                     data=scores_basal,
                     family=binomial)
    coef <- summary(lmres)$coefficients$cond
    model_p[i,"Score_Logit"] <- coef[i,"Estimate"]
    model_p[i,"Score_P"] <- coef[i,"Pr(>|z|)"]
    model_p[i,"Smoke_Logit"] <- coef["Smoking_StatusCurrent","Estimate"]
    model_p[i,"Smoke_P"] <- coef["Smoking_StatusCurrent","Pr(>|z|)"]
  }
  model_p$Score_Q <- p.adjust(model_p$Score_P, method = "fdr")
  model_p$Smoke_Q <- p.adjust(model_p$Smoke_P, method = "fdr")
  list(model_results = model_p, module_scores = scores_basal)
}
glmm_module_res_lgbvshgb <- vam_model_logistic(scores = vam_scores)
lgbvshgb_model_results <- glmm_module_res_lgbvshgb$model_results %>% arrange(Score_Q)

## Figure 2G: Violins for Select Outside Modules in LGB vs. HGB Cells
sig_diff_interesting <- c("Mascaux_M1", "Mascaux_M2", "Merrick_HG_Up", "Merrick_HG_Down", "bx_lightgreen_pos", "bx_blue_pos")
celltype_cols <- readRDS("celltype_colors.rds")
sig_diff_violin <- lapply(sig_diff_interesting, function(x) {
  plt <- ggplot(scores_basal, aes_string(x = "KCluster", y = x, fill = "CellType")) + 
    geom_violin(scale = "width") + xlab("Cluster") + ylab("VAM Score") + ggtitle(x) + 
    theme_classic() + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + scale_fill_manual(values = celltype_cols) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")
})
sig_diff_interesting_violin <- plot_grid(plotlist = sig_diff_violin, align = "hv", ncol = 2)
pdf("LGB_HGB_OutsideModules_Violins.pdf", width = 6, height = 6)
sig_diff_interesting_violin
dev.off()

## Compare relative proportions of Basal Sample Grade (High- vs. Low-Grade Samples) and Smoking Status (Current vs. Smoker) in non-basal cell types
sce_hasenoughcells <- subsetSCECols(sce, index = which(sce$BasalGroupMixed %in% c("Low Grade Basal Sample", "High Grade Basal Sample")))
sce_hasenoughcells[["BasalGroupMixed"]] <- factor(sce_hasenoughcells[["BasalGroupMixed"]], levels = c("Low Grade Basal Sample", "High Grade Basal Sample"), ordered = FALSE)
colData(sce_hasenoughcells) <- within(colData(sce_hasenoughcells), BasalGroupMixed <- relevel(BasalGroupMixed, ref = "Low Grade Basal Sample"))
sce_hasenoughcells[["Smoking_Status"]] <- factor(sce_hasenoughcells[["Smoking_Status"]], levels = c("Former", "Current"), ordered = FALSE)
colData(sce_hasenoughcells) <- within(colData(sce_hasenoughcells), Smoking_Status <- relevel(Smoking_Status, ref = "Former"))
notbasalcells <- which(!sce_hasenoughcells$CellType %in% c("Basal Cells", "KRT5+/MUC5B+ Cells", "KRT5+/SCGB1A1+ Cells"))
sce_nobasalcomp <- subsetSCECols(sce_hasenoughcells, index = notbasalcells)
sce_nobasalcomp$KCluster <- factor(sce_nobasalcomp$KCluster,
                                   levels = paste0("K",12:35),
                                   ordered = FALSE)

mixed_cellshifts_enoughcells_percluster_wsite <- function(sce_object) {
  model_with_factor_association <- 
    sce_object |>
    sccomp_estimate( 
      formula_composition = ~ 1 + BasalGroupMixed + Smoking_Status + PCGA02_site, 
      .sample = Sample, 
      .cell_group = KCluster, 
      inference_method = "hmc",
      enable_loo = TRUE
    ) |> 
    sccomp_remove_outliers(cores = 4) |> 
    sccomp_test()
  plot_sum_basalgroup_smoking <- model_with_factor_association |> 
    sccomp_boxplot(factor = "BasalGroupMixed")
  list(fullmodel = model_with_factor_association, fullmodel_plot = plot_sum_basalgroup_smoking)
}
nonimmune_basalgroup_percluster_wsite <- mixed_cellshifts_enoughcells_percluster_wsite(sce_object = sce_nobasalcomp)
nonimmune_basalgroup_model_res <- nonimmune_basalgroup_percluster_wsite$fullmodel %>% arrange(parameter, KCluster) %>% select(KCluster:c_FDR)

## Figure 2H: Cell Type x Basal Sample Grade and Smoking Status model results bubble plot
stats_values <- nonimmune_basalgroup_model_res %>% filter(factor %in% c("BasalGroupMixed", "Smoking_Status"))
stats_values$c_FDR <- as.numeric(as.vector(stats_values$c_FDR))
stats_values$c_FDR[stats_values$c_FDR == 0] <- 2e-16
stats_values$c_effect[which(stats_values$c_effect < -3)] <- -3
stats_values$c_effect[which(stats_values$c_effect > 3)] <- 3
stats_values$logp <- -log10(stats_values$c_FDR)
stats_values$sig <- NA
for (i in 1:length(stats_values$c_FDR)) {
  if (is.na(stats_values$c_FDR[i])) 
    stats_values$sig[i] <- "" 
  else if (stats_values$c_FDR[i] > 0.05)
    stats_values$sig[i] <- ""
  else if (stats_values$c_FDR[i] > 0.01) 
    stats_values$sig[i] <- "*" 
  else if (stats_values$c_FDR[i] > 0.001) 
    stats_values$sig[i] <- "**"  
  else stats_values$sig[i] <- "***"
}

stats_values$parameter <- gsub("BasalGroupMixedHigh Grade Basal Sample", "High Grade vs. Low Grade Basal Sample", stats_values$parameter)
stats_values$parameter <- gsub("Smoking_StatusCurrent", "Current vs. Former Smoker", stats_values$parameter)
stats_values$parameter <- factor(stats_values$parameter,
                                 levels = c("High Grade vs. Low Grade Basal Sample", "Current vs. Former Smoker"),
                                 ordered = TRUE)

colors <- colorRampPalette(c("#053061", "#053061", "#053061", "#2166AC", "#4393C3", 
                             "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#67001F", 
                             "#67001F", "#67001F"))(200)
cellshift_cluster_basalgroup_bubble <- ggplot(stats_values, aes(x = parameter, y = KCluster, label = sig, color = c_effect, 
                                                                size = ifelse(is.na(logp), 2, logp), shape = ifelse(is.na(c_effect), "Missing", "Present"))) + 
  geom_point(pch=1, colour = "black", stroke = 2) +
  geom_point(stat = "identity", stroke = 2, na.rm = T, shape = 19) +  
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "", y = "", color = "Effect Size", size = "-log10(FDR)") + 
  guides(size = "legend", color = "colorbar") +
  geom_text(color = "white", size = 3, vjust = 0.9) +
  scale_x_discrete(labels = c("High Grade vs. Low Grade\nBasal Sample", "Current vs.\nFormer Smoker")) +
  scale_y_discrete(limits = rev) +
  scale_shape_manual(values = c(Missing = 4, Present = 19)) +
  scale_color_gradientn(colours = colors, limits = c(-3, 3), breaks = c(-3, 0, 3))

celltype_colorbar_df <- data.frame(table(sce$KCluster, sce$CellType))
celltype_colorbar_df <- celltype_colorbar_df[which(celltype_colorbar_df$Freq > 0),]
celltype_colorbar_df$X <- 1
celltype_colorbar_df$Freq <- c()
colnames(celltype_colorbar_df) <- c("KCluster", "CellType", "X")
celltype_colorbar_df <- celltype_colorbar_df[which(!celltype_colorbar_df$KCluster %in% paste0("K",1:11)),]
celltype_colorbar_df$KCluster <- factor(celltype_colorbar_df$KCluster,
                                        levels = paste0("K",35:12),
                                        ordered = TRUE)
celltype_colorbar_ggbar <- ggplot(celltype_colorbar_df, aes(x = X, y = KCluster, fill = CellType)) + geom_bar(position="fill", stat="identity") +
  theme_minimal() + scale_fill_manual(values = celltype_colors) + xlab("") + ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_blank(), legend.position = "none")

df_basalgroup <- data.frame(table(sce$KCluster, sce$BasalGroupMixed))
df_basalgroup <- df_basalgroup[which(!df_basalgroup$Var1 %in% paste0("K",1:11)),]
df_basalgroup$Var1 <- factor(df_basalgroup$Var1, levels = paste0("K",35:12), ordered = TRUE)
df_basalgroup$Var2 <- factor(df_basalgroup$Var2, levels = rev(c("Low Grade Basal Sample", "High Grade Basal Sample")), ordered = TRUE)
basalgroup_colors <- c("Low Grade Basal Sample" = "blue",
                       "High Grade Basal Sample" = "red")
basalgroup_sbp <- ggplot(df_basalgroup, aes(fill=Var2, x=Freq, y=Var1)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() + 
  ylab("") + xlab("% Basal Sample") + ggtitle("") +
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = scales::percent, expand = expansion(add = c(0,0.02))) + 
  scale_fill_manual(limits = rev, name = "Sample Type", values = basalgroup_colors) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")) 

df_smoke <- data.frame(table(sce$KCluster, sce$Smoking_Status))
df_smoke <- df_smoke[which(!df_smoke$Var1 %in% paste0("K",1:11)),]
df_smoke$Var1 <- factor(df_smoke$Var1, levels = paste0("K",35:12), ordered = TRUE)
df_smoke$Var2 <- factor(df_smoke$Var2, levels = rev(levels(sce$Smoking_Status)), ordered = TRUE)
smoking_sbp <- ggplot(df_smoke, aes(fill=Var2, x=Freq, y=Var1)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() + 
  ylab("") + xlab("% Smoking") + ggtitle("") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = scales::percent, expand = expansion(add = c(0,0.02))) + scale_fill_manual(name = "Smoking Status", values = smoke_cols) +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")) 

pdf("SCComp_Nonimmune_BubblePlot_PerCluster_wSite.pdf", width = 12, height = 6)
plot_grid(plotlist = list(celltype_colorbar_ggbar, basalgroup_sbp, smoking_sbp, cellshift_cluster_basalgroup_bubble), align = "h", nrow = 1, rel_widths = c(1, 7, 7, 7))
dev.off()

## Figure 2I: Select significantly different cluster proportions by Basal Sample Grade
notbasalcells <- which(!sce_hasenoughcells$CellType %in% c("Basal Cells", "KRT5+/MUC5B+ Cells", "KRT5+/SCGB1A1+ Cells"))
sce_nobasalcomp <- subsetSCECols(sce_hasenoughcells, index = notbasalcells)
sce_nobasalcomp$KCluster <- factor(sce_nobasalcomp$KCluster,
                                   levels = paste0("K",12:35),
                                   ordered = FALSE)
percent_celltype <- data.frame(prop.table(table(sce_nobasalcomp$Sample, sce_nobasalcomp$KCluster), margin = 1))
colnames(percent_celltype) <- c("Sample", "KCluster", "Proportion")
sample_metadata_df <- data.frame(unique(data.frame("Sample" = sce_nobasalcomp$Sample, "SmokingStatus" = sce_nobasalcomp$Smoking_Status, "BasalGroup" = sce_nobasalcomp$BasalGroupMixed)))
percent_celltype <- left_join(percent_celltype, 
                              sample_metadata_df,
                              by = c("Sample"))
cluster_metadata_df <- data.frame(unique(data.frame("KCluster" = sce$KCluster, "CellType" = sce$CellType)))
cluster_metadata_df$KCluster <- factor(cluster_metadata_df$KCluster, levels = levels(cluster_metadata_df$KCluster), ordered = FALSE)
percent_celltype <- left_join(percent_celltype, 
                              cluster_metadata_df,
                              by = c("KCluster"))
percent_celltype$KCluster_CellType <- paste(percent_celltype$KCluster, percent_celltype$CellType, sep = " - ")
percent_celltype$BasalGroup <- factor(percent_celltype$BasalGroup,
                                      levels = c("Low Grade Basal Sample", "High Grade Basal Sample"),
                                      ordered = TRUE)

cellshifts_basalgroup_fullmodel_mixed_violin_interesting_bycluster <- lapply(c("K21 - Club Cells", "K25 - Ciliated Cells",
                                                                               "K31 - Airway Smooth Muscle Cells", "K32 - Fibroblasts"), FUN = function(x) {
                                                                                 df <- percent_celltype[which(percent_celltype$KCluster_CellType == x),]
                                                                                 p <- ggplot(df, aes(x = BasalGroup, y = Proportion, fill = SmokingStatus)) + geom_boxplot() +
                                                                                   theme_classic() + 
                                                                                   scale_x_discrete(labels = c("Low Grade\nBasal Sample", "High Grade\nBasal Sample")) +
                                                                                   scale_y_continuous(labels = scales::percent) +
                                                                                   scale_fill_manual(values = smoke_cols) +
                                                                                   labs(title = x, y = "% Non-Basal Cells", x = "", fill = "Smoking Status") +
                                                                                   theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
                                                                                 p
                                                                               })
pdf("cellshifts_basalgroup_nomixed_violin_significant_bycluster_wsite.pdf", width = 8, height = 8)
plot_grid(plotlist = cellshifts_basalgroup_fullmodel_mixed_violin_interesting_bycluster, ncol = 2, align = "h")
dev.off()

########## ------------- SUPPLEMENTARY FIGURES/TABLES, FIGURE 2 ------------- ##########

## Supplementary Figure 4: Nonimmune Biopsy Cell Cluster UMAP
umap_coords <- reducedDim(altExp(sce), "celda_UMAP")
colnames(umap_coords) <- c("UMAP1", "UMAP2")
cluster_df <- cbind(data.frame(Clusters = paste0("K", celdaClusters(sce))), umap_coords)
cluster_df$Clusters <- factor(cluster_df$Clusters, levels = paste0("K",levels(celdaClusters(sce))), ordered = TRUE)
cluster_celltype_df <- cbind(data.frame(Clusters = paste0(sce$CellType, " (K", celdaClusters(sce),")")), umap_coords)
celltype_df <- cbind(data.frame(CellType = sce$CellType), umap_coords)
cluster_colors <- distinctColors(
  n = length(unique(celdaClusters(sce))),
  hues = c("red", "cyan", "orange", "pink1", "blue", "yellow", "purple", "green", "forestgreen", "magenta"))

cluster_umap <- ggplot(cluster_df, aes(x = UMAP1, y = UMAP2, color = Clusters)) + geom_point(size = 0.25) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                          axis.ticks = element_blank()) + scale_color_manual(values = cluster_colors)
add_labels <- function(df, label_column, xcoord = "UMAP1", ycoord = "UMAP2", plt) {
  centroidList <- lapply(unique(df[,label_column]), function(x) {
    df.sub <- df[df[,label_column] == x, ]
    median1 <- stats::median(df.sub[, xcoord])
    median2 <- stats::median(df.sub[, ycoord])
    data.frame(median1 = median1,
               median2 = median2,
               x = x)
  })
  centroid <- do.call(rbind, centroidList)
  centroid <- data.frame(
    Dimension_1 = as.numeric(centroid[, 1]),
    Dimension_2 = as.numeric(centroid[, 2]),
    Clusters = centroid[, 3]
  )
  colnames(centroid)[seq(2)] <- c(xcoord, ycoord)
  
  label_object <- plt + geom_point(
    data = centroid,
    mapping = aes_string(x = xcoord,y = ycoord),
    size = 0,
    alpha = 0,
    inherit.aes = FALSE
  ) +
    ggrepel::geom_text_repel(
      data = centroid,
      mapping = ggplot2::aes_string(label = label_column),
      size = 3.5,
      max.overlaps = Inf,
      force = 5,
      force_pull = 0.5,
      show.legend = FALSE,
      color = "black"
    )
  label_object
}
centroid_umap_labeled <- add_labels(cluster_df, "Clusters", plt = cluster_umap)

alphabetical_cluster_celltype_levels <- unique(cluster_celltype_df$Clusters)
cluster_celltype_umap_levels <- c()
for(i in levels(cluster_df$Clusters)) {
  index <- grep(paste0(i,")"), alphabetical_cluster_celltype_levels)
  cluster_celltype_umap_levels <- c(cluster_celltype_umap_levels, alphabetical_cluster_celltype_levels[index])
}
cluster_celltype_df$Clusters <- factor(cluster_celltype_df$Clusters, levels = cluster_celltype_umap_levels, ordered = TRUE)
cluster_celltype_umap <- ggplot(cluster_celltype_df, aes(x = UMAP1, y = UMAP2, color = Clusters)) + geom_point(size = 0.25) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                          axis.ticks = element_blank()) + scale_color_manual(values = cluster_colors)
cluster_celltype_umap_labeled <- add_labels(cluster_celltype_df, "Clusters", plt = cluster_celltype_umap)
cluster_celltype_umap_labeled <- cluster_celltype_umap_labeled + labs(color = "Cluster")
pdf("NonimmuneCell_CellType+Cluster_UMAP.pdf", width = 8, height = 6)
cluster_celltype_umap_labeled
dev.off()

## Supplementary Figure 5: Nonimmune Probability Map
mod_map <- celdaProbabilityMap(sce, useAssay = "decontXcounts")
pdf("NonimmuneCell_ProbabilityMap.pdf", width = 10, height = 10)
mod_map
dev.off()

## Supplementary Figure 6: Low Grade Heterogeneity Modules
modules <- read.table("070924_Nonimmune_celda_K35_L73_CeldaEdit.tsv", sep = "\t", header = TRUE)
modules <- lapply(as.list(modules), function(x) {
  x <- x[x != ""]
  x
})
modules_list <- as.list(modules)
marker_mods <- c("L22", "L70", "L72", "L60", "L58", "L59", "L48", "L49", "L51", "L26", "L27", "L28",
                 "L52", "L53", "L54", "L55", "L30", "L62", "L63")
markers <- c("KRT5", "CEACAM5", "KRT8", "MUC5AC", "SCGB1A1", "SCGB3A1", "BPIFA1", "BPIFB1", "FOXJ1",
             "TUBA1B", "GSTP1", "ANXA2", "SNTN", "TPPP3/TUBA1A", "DYNLL1", "GSTA1-2",
             "AKR1C1-3", "ALDH3A1/ADH7", "CYP1A1/CYP1B1")
m <- match(marker_mods, rownames(sce_modular))
scaled_expr <- t(assay(sce_modular,"LogNormalize")[m,])
scaled_expr <- scale(scaled_expr)
colnames(scaled_expr) <- paste0(marker_mods," (",markers,")")
df <- cbind(data.frame(Cluster = sce$KCluster), scaled_expr)
df_markersummary <- df %>% group_by(Cluster) %>% summarize_all(mean)
df_markersummary_pivot <- df_markersummary %>% pivot_longer(!Cluster, names_to = "Module", values_to = "Expression")
df_markersummary_pivot$Cluster <- factor(df_markersummary_pivot$Cluster,
                                         levels = rev(paste0("K",1:35)),
                                         ordered = TRUE)
df_markersummary_pivot$Module <- factor(df_markersummary_pivot$Module, levels = colnames(scaled_expr), ordered = TRUE)
midpoint <- mean(df_markersummary_pivot$Expression)
marker_scalehm <- ggplot(df_markersummary_pivot, aes(x = Module, y = Cluster, fill = Expression)) + geom_tile() +
  theme_minimal() + ylab("") + xlab("") +
  scale_fill_gradient2(name = "Average\nScaled\nNormalized\nExpression", midpoint = midpoint,
                       low = "blue", mid = "white", high = "red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, -0.5), "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.ticks = element_blank())
marker_scalehm

celltype_colorbar_df <- data.frame(Cluster = levels(sce$CellType_Cluster), X = 1)
celltype_colorbar_df$Cluster <- factor(celltype_colorbar_df$Cluster,
                                       levels = rev(levels(sce$CellType_Cluster)), ordered = TRUE)
celltype_cluster_cols <- c(rep(celltype_colors[["Basal Cells"]], 9), celltype_colors[["KRT5+/MUC5B+ Cells"]], celltype_colors[["KRT5+/SCGB1A1+ Cells"]],
                           celltype_colors[["Perigoblet Cells"]], rep(celltype_colors[["Goblet Cells"]], 5), rep(celltype_colors[["Perigoblet Cells"]], 5),
                           celltype_colors[["Mucous SMG Cells"]], celltype_colors[["Serous SMG Cells"]], rep(celltype_colors[["Perigoblet Cells"]], 6),
                           celltype_colors[["Airway Smooth Muscle Cells"]], rep(celltype_colors[["Fibroblasts"]], 2), rep(celltype_colors[["Endothelial Cells"]], 2))
names(celltype_cluster_cols) <- rev(levels(celltype_colorbar_df$Cluster))
celltype_colorbar_ggbar <- ggplot(celltype_colorbar_df, aes(x = X, y = Cluster, fill = Cluster)) + geom_bar(position="fill", stat="identity") +
  theme_minimal() + scale_fill_manual(values = celltype_cluster_cols) + xlab("") + ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_blank(), legend.position = "none")

df_hist <- data.frame(table(sce$KCluster, sce$Histology))
df_hist$Var1 <- factor(df_hist$Var1, levels = paste0("K",35:1), ordered = TRUE)
df_hist$Var2 <- factor(df_hist$Var2, levels = rev(levels(sce$Histology)), ordered = TRUE)
histology_sbp <- ggplot(df_hist, aes(fill=Var2, x=Freq, y=Var1)) +
  geom_bar(position="fill", stat="identity") + theme_classic() +
  ylab("") + xlab("% Histology") + ggtitle("") +
  scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), expand = expansion(add = c(0,0.02))) + scale_fill_manual(name = "Histology", values = hist_cols) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_legend(nrow = 5, reverse = TRUE))

df_smoke <- data.frame(table(sce$KCluster, sce$Smoking_Status))
df_smoke$Var1 <- factor(df_smoke$Var1, levels = paste0("K",35:1), ordered = TRUE)
df_smoke$Var2 <- factor(df_smoke$Var2, levels = rev(levels(sce$Smoking_Status)), ordered = TRUE)
smoke_sbp <- ggplot(df_smoke, aes(fill=Var2, x=Freq, y=Var1)) +
  geom_bar(position="fill", stat="identity") + theme_classic() +
  ylab("") + xlab("% Smoking") + ggtitle("") +
  scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), expand = expansion(add = c(0,0.02))) + scale_fill_manual(name = "Smoking Status", values = smoke_cols) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_legend(nrow = 2, reverse = TRUE))

stats_values <- nonimmune_model_res %>% filter(factor %in% c("HistologyBinary", "Smoking_Status"))
stats_values$c_FDR[stats_values$c_FDR == 0] <- 2e-16
stats_values$c_effect[which(stats_values$c_effect < -3)] <- -3 
stats_values$c_effect[which(stats_values$c_effect > 3)] <- 3 
stats_values$logp <- -log10(stats_values$c_FDR)
stats_values$sig <- NA
for (i in 1:length(stats_values$c_FDR)) {
  if (is.na(stats_values$c_FDR[i])) 
    stats_values$sig[i] <- "" 
  else if (stats_values$c_FDR[i] > 0.05)
    stats_values$sig[i] <- ""
  else if (stats_values$c_FDR[i] > 0.01) 
    stats_values$sig[i] <- "*" 
  else if (stats_values$c_FDR[i] > 0.001) 
    stats_values$sig[i] <- "**"  
  else stats_values$sig[i] <- "***"
}

stats_values$parameter <- gsub("HistologyBinaryHigh Grade", "Severe Dysplasia - LUSC vs. Normal - Moderate Dysplasia", stats_values$parameter)
stats_values$parameter <- gsub("Smoking_StatusCurrent", "Current vs. Former Smoker", stats_values$parameter)
stats_values$parameter <- factor(stats_values$parameter,
                                 levels = c("Severe Dysplasia - LUSC vs. Normal - Moderate Dysplasia", "Current vs. Former Smoker"),
                                 ordered = TRUE)
stats_values$KCluster <- factor(stats_values$KCluster,
                                levels = paste0("K",35:1), # Change to 11 if just plotting basal cells as below.
                                ordered = TRUE)

colors <- colorRampPalette(c("#053061", "#053061", "#053061", "#2166AC", "#4393C3", 
                             "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#67001F", 
                             "#67001F", "#67001F"))(200)
cellshift_cluster_histology_bubble <- ggplot(stats_values, aes(x = parameter, y = KCluster, label = sig, color = c_effect, 
                                                               size = ifelse(is.na(logp), 2, logp), shape = ifelse(is.na(c_effect), "Missing", "Present"))) + 
  geom_point(pch=1, colour = "black", stroke = 2) +
  geom_point(stat = "identity", stroke = 2, na.rm = T, shape = 19) +  
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5),
        axis.ticks = element_blank()) +
  labs(x = "", y = "", color = "Effect Size", size = "-log10(FDR)") + 
  guides(size = "legend", color = "colorbar") +
  geom_text(color = "white", size = 3, vjust = 0.9) +
  scale_x_discrete(labels = c("Severe Dysplasia - LUSC\nvs.\nNormal - Moderate Dysplasia", "Current vs.\nFormer Smoker")) +
  scale_shape_manual(values = c(Missing = 4, Present = 19)) +
  scale_color_gradientn(colours = colors, limits = c(-3, 3), breaks = c(-3, 0, 3)) 

# Combine HM, boxplots, and sccomp figure
pdf("LowGrade_Histology_Associated_Modules_HM_Nonimmune_wBoxPlot+SCComp.pdf", width = 15, height = 8)
plot_grid(celltype_colorbar_ggbar, marker_scalehm, histology_sbp, smoke_sbp, cellshift_cluster_histology_bubble, align = "h", nrow = 1, rel_widths = c(2,5,2,2,2))
dev.off()

## Supplementary Figure 7: LGB vs. HGB Module Volcano Plot
volcano_df$Label <- ifelse((abs(volcano_df$CellType_Beta) > 1) & (volcano_df$CellType_Q < 0.05), "Label", "NoLabel")
label_cols <- c("Label" = "red", "NoLabel" = "black")
volcano_gg <- ggplot(volcano_df, aes(x = CellType_Beta, y = -log10(CellType_Q), color = Label)) + geom_point() + 
  geom_text_repel(aes(label=ifelse(Label == "Label",as.character(Module),'')), show.legend = FALSE) + 
  theme_classic() + xlim(-4, 4) + scale_color_manual(values = label_cols) + 
  guides(color = "none") + 
  xlab("Beta, High-Grade - Low-Grade Basal") + ylab("-log10(FDR)") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "blue") + geom_vline(xintercept=-1, linetype="dashed", color = "blue") + 
  geom_vline(xintercept=1, linetype="dashed", color = "blue") 
pdf("LGB_HGB_Module_Volcano_wMixed_wSite.pdf", width = 4, height = 4)
volcano_gg
dev.off()

## Supplementary Figure 8: Outside Modules in LGB vs. HGB Volcano Plot
volcano_df <- read.xlsx("LGB_HGB_OutsideModules_VAM_LogR_wSite_FullModel.xlsx")
label_cols <- c("Label" = "red", "NoLabel" = "black")
volcano_gg <- ggplot(volcano_df, aes(x = Score_Logit, y = -log10(Score_Q), color = Label)) + geom_point() + 
  geom_text_repel(aes(label=ifelse(Label == "Label",as.character(Gene.Set),'')), show.legend = FALSE) + 
  theme_classic() + xlim(-10, 10) + scale_color_manual(values = label_cols) + 
  guides(color = "none") + #theme(legend.position = "none") +
  xlab("Beta, High-Grade - Low-Grade Basal") + ylab("-log10(FDR)") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "blue") + geom_vline(xintercept=-3, linetype="dashed", color = "blue") + 
  geom_vline(xintercept=3, linetype="dashed", color = "blue") 
pdf("LGB_HGB_OutsideModule_Volcano_wSite.pdf", width = 4, height = 4)
volcano_gg
dev.off()

## Supplementary Table 7: Enrichment Terms for Nonimmune Modules (EnrichR)
modules <- read.table("070924_Nonimmune_celda_K35_L73_CeldaEdit.tsv", sep = "\t", header = TRUE)
modules <- lapply(as.list(modules), function(x) {
  x <- x[x != ""]
  x
})
modules_enrichr <- lapply(modules, function(x) {
  enrichments <- enrichr(genes = x, databases = "GO_Biological_Process_2023")
  enrichments
})
modules_enrichr_df <- data.frame(Module = paste0("L",1:73), 'Top 10 Enriched Genes' = rep("",73), 'Top 5 Enrichment Terms' = rep("",73), check.names = FALSE)
for(i in seq_along(modules_enrichr_df$Module)) {
  df <- modules_enrichr[[i]][["GO_Biological_Process_2023"]]
  genes <- paste(modules_enrichr[[i]][["GO_Biological_Process_2023"]][["Genes"]][1:5], collapse = ";")
  genes <- strsplit(genes, split = ";")[[1]]
  genes <- unique(genes)
  genes <- genes[1:min(length(genes),10)]
  modules_enrichr_df$`Top 10 Enriched Genes`[i] <- paste(genes, collapse = ";")
  modules_enrichr_df$`Top 5 Enrichment Terms`[i] <- paste(modules_enrichr[[i]][["GO_Biological_Process_2023"]][["Term"]][1:5], collapse = ";")
}
write.xlsx(modules_enrichr_df, file = "Nonimmune_Module_Enrichr.xlsx")

## Supplementary Table 8: Cell Type x Histology and Smoking Status (SCComp)
write.xlsx(nonimmune_model_res, file = "NonimmuneClusterbyHistology_wSiteEffect_FullModel.xlsx", rowNames = FALSE) 

## Supplementary Table 9: LGB vs. HGB Module Results
write.xlsx(volcano_df, "LGB_vs_HGB_Module_Analysis_wSite_031225.xlsx")

## Supplementary Table 10: Outside Modules in LGB vs. HGB Module Results
write.xlsx(lgbvshgb_model_results, "LGB_HGB_OutsideModules_VAM_LogR_wSite_FullModel.xlsx", rowNames = FALSE)

## Supplementary Table 11: Cell Type Clusters x Basal Sample Grade and Smoking Status (SCComp)
write.xlsx(nonimmune_basalgroup_model_res, file = "NonBasalClusterbyBasalGroup_wSiteEffect_FullModel.xlsx", rowNames = FALSE) 

saveRDS(sce_full, "AllBiopsyCell_SCE_Final.rds")
saveRDS(sce, "NonimmuneBiopsyCell_SCE_Final")
