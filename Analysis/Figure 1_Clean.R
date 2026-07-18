library(celda)
library(singleCellTK)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(RColorBrewer)
library(sccomp)
library(cmdstanr)
library(stringr)
library(magrittr)
library(loo)
library(openxlsx)

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/")

# I updated some metadata and added BasalGroup from Figure2_Clean.R
sce <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/AllBiopsyCell_SCE_Final.rds")

sce$HistologyFull <- sce$Histology
sce$Histology <- as.character(sce$Histology)
sce$Histology[which(sce$Histology == "Denuded Bronchial Mucosa")] <- "Normal"
sce$Histology[which(sce$Histology == "Inflammation")] <- "Normal"
sce$Histology <- factor(sce$Histology, 
                        levels = c("Normal", "Basal Cell Hyperplasia", "Metaplasia", "Mild Dysplasia",
                                   "Moderate Dysplasia", "Severe Dysplasia", "CIS", "LUSC"),
                        ordered = TRUE)

celltype_colors <- readRDS("celltype_colors.rds") 
hist_cols <- readRDS("histology_colors.rds")
smoke_cols <- readRDS("smoking_colors.rds")

f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")), 
                                    colData = colData(altExp(sce)))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")

## Figure 1B: Cell Type UMAP
sce$CellType_Cluster <- paste0(sce$CellType, "_K", sce$celda_cell_cluster)
sce$CellType_ClusterRange <- as.character(sce$CellType_Cluster)
k1.3 <- which(sce$celda_cell_cluster %in% c(1:3)); sce$CellType_ClusterRange[k1.3] <- "KRT5+ Basal Cells (K1-K3)"
k4.7 <- which(sce$celda_cell_cluster %in% c(4:7)); sce$CellType_ClusterRange[k4.7] <- "KRT5+ Basal Cells (K4-K7)"
k8.9 <- which(sce$celda_cell_cluster %in% c(8:9)); sce$CellType_ClusterRange[k8.9] <- "MUC5AC+ Goblet Cells (K8-K9)"
k10.11 <- which(sce$celda_cell_cluster %in% c(10:11)); sce$CellType_ClusterRange[k10.11] <- "SCGB1A1+ Club Cells (K10-K11)"
k12 <- which(sce$celda_cell_cluster %in% c(12)); sce$CellType_ClusterRange[k12] <- "MUC5B+/PRB1+ SMG Cells (K12)"
k13.15 <- which(sce$celda_cell_cluster %in% c(13:15)); sce$CellType_ClusterRange[k13.15] <- "FOXJ1+ Ciliated Cells (K13-K15)"
k16.17 <- which(sce$celda_cell_cluster %in% c(16:17)); sce$CellType_ClusterRange[k16.17] <- "COL1A1+ Fibroblasts (K16-K17)"
k18 <- which(sce$celda_cell_cluster %in% c(18)); sce$CellType_ClusterRange[k18] <- "VWF+ Endothelial Cells (K18)"
k19.21 <- which(sce$celda_cell_cluster %in% c(19:21)); sce$CellType_ClusterRange[k19.21] <- "CD8+ T Cells (K19-K21)"
k22.23 <- which(sce$celda_cell_cluster %in% c(22:23)); sce$CellType_ClusterRange[k22.23] <- "CD4+ T Cells (K22-K23)"
k24 <- which(sce$celda_cell_cluster %in% c(24)); sce$CellType_ClusterRange[k24] <- "MS4A1+ B Cells (K24)"
k25 <- which(sce$celda_cell_cluster %in% c(25)); sce$CellType_ClusterRange[k25] <- "MZB1+ Plasma Cells (K25)"
k26.27 <- which(sce$celda_cell_cluster %in% c(26:27)); sce$CellType_ClusterRange[k26.27] <- "CD68+ Macrophages (K26-K27)"
k28 <- which(sce$celda_cell_cluster %in% c(28)); sce$CellType_ClusterRange[k28] <- "FCER1A+ Dendritic Cells (K28)"
k29 <- which(sce$celda_cell_cluster %in% c(29)); sce$CellType_ClusterRange[k29] <- "CSF3R+ Neutrophils (K29)"
k30 <- which(sce$celda_cell_cluster %in% c(30)); sce$CellType_ClusterRange[k30] <- "MS4A2+ Mast Cells (K30)"
altExp(sce)$CellType_ClusterRange <- sce$CellType_ClusterRange

celltype_rangecols <- rep(NA, length(unique(sce$CellType_ClusterRange))); names(celltype_rangecols) <- unique(sce$CellType_ClusterRange)
celltype_rangecols["KRT5+ Basal Cells (K1-K3)"] <- celltype_colors[["Basal Cells"]]; celltype_rangecols["KRT5+ Basal Cells (K4-K7)"] <- celltype_colors[["Basal Cells"]]; 
celltype_rangecols["MUC5AC+ Goblet Cells (K8-K9)"] <- celltype_colors[["Goblet Cells"]]; celltype_rangecols["SCGB1A1+ Club Cells (K10-K11)"] <- celltype_colors[["Club Cells"]];  
celltype_rangecols["MUC5B+/PRB1+ SMG Cells (K12)"] <- celltype_colors[["SMG Cells"]]; celltype_rangecols["FOXJ1+ Ciliated Cells (K13-K15)"] <- celltype_colors[["Ciliated Cells"]]; 
celltype_rangecols["COL1A1+ Fibroblasts (K16-K17)"] <- celltype_colors[["Fibroblasts"]]; celltype_rangecols["VWF+ Endothelial Cells (K18)"] <- celltype_colors[["Endothelial Cells"]]; 
celltype_rangecols["CD8+ T Cells (K19-K21)"] <- celltype_colors[["CD8+ T Cells"]]; celltype_rangecols["CD4+ T Cells (K22-K23)"] <- celltype_colors[["CD4+ T Cells"]]; 
celltype_rangecols["MS4A1+ B Cells (K24)"] <- celltype_colors[["B Cells"]]; celltype_rangecols["MZB1+ Plasma (K25)"] <- celltype_colors[["Plasma Cells"]];  
celltype_rangecols["CD68+ Macrophages (K26-K27)"] <- celltype_colors[["Macrophages"]];  celltype_rangecols["FCER1A+ Dendritic Cells (K28)"] <- celltype_colors[["Dendritic Cells"]];  
celltype_rangecols["CSF3R+ Neutrophils (K29)"] <- celltype_colors[["Neutrophils"]];  celltype_rangecols["MS4A2+ Mast Cells (K30)"] <- celltype_colors[["Mast Cells"]]

celltype_umap <- plotSCEDimReduceColData(altExp(sce), colorBy = "CellType_ClusterRange", reducedDimName = "celda_UMAP", dotSize = 0.25, xlab = "UMAP1", ylab = "UMAP2")
celltype_umap_altered <- celltype_umap + theme_classic() + scale_color_manual(values = celltype_rangecols) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), legend.position = "none")
pdf("AllCell_CellType_UMAP_wMarkers.pdf", width = 3, height = 3)
celltype_umap_altered
dev.off()

smoking_df <- cbind(data.frame(Smoking = sce$Smoking_Status), umap_coords)
smoking_umap <- ggplot(smoking_df, aes(x = UMAP1, y = UMAP2, color = Smoking)) + geom_point(size = 0.25) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                          axis.ticks = element_blank()) + scale_color_manual(name = "Smoking Status", values = smoke_cols, na.value = "gray80")
pdf("All_Cell_Smoking_UMAP.pdf", width = 4.5, height = 3)
smoking_umap
dev.off()

histology_df <- cbind(data.frame(Histology = sce$Histology), umap_coords)
hist_umap <- ggplot(histology_df, aes(x = UMAP1, y = UMAP2, color = Histology)) + geom_point(size = 0.25) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                          axis.ticks = element_blank()) + scale_color_manual(values = hist_cols, na.value = "gray80")
pdf("All_Cell_CollapsedHistology_UMAP.pdf", width = 5, height = 3)
hist_umap
dev.off()

## Figure 1C: Marker Module HM
markers <- c("KRT5", "MUC5AC", "SCGB1A1", "MUC5B", "PRB1", "FOXJ1",  "COL1A1", "VWF", 
             "CD8A", "CD4", "MS4A1", "MZB1", "CD68", "FCER1A", "CSF3R", 
             "MS4A2")
modules <- read.table("070124_celda_K30_L104.tsv", sep = "\t", header = TRUE)
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
colnames(scaled_expr)[4] <- "L90 (MUC5B/PRB1)"
scaled_expr <- scaled_expr[,c(setdiff(1:16,c(5)))]
sce_modular$CellType <- sce$CellType
df <- cbind(data.frame(CellType = sce_modular$CellType), scaled_expr)
df_markersummary <- df %>% group_by(CellType) %>% summarize_all(mean)
df_markersummary_pivot <- df_markersummary %>% pivot_longer(!CellType, names_to = "Module", values_to = "Expression")
df_markersummary_pivot$CellType <- factor(df_markersummary_pivot$CellType,
                                          levels = c("Basal Cells", "Goblet Cells", "Club Cells",
                                                     "SMG Cells", "Ciliated Cells", "Fibroblasts", "Endothelial Cells", 
                                                     "CD8+ T Cells", "CD4+ T Cells", "B Cells", "Plasma Cells", "Macrophages", 
                                                     "Dendritic Cells", "Neutrophils", "Mast Cells"),
                                          ordered = TRUE)
df_markersummary_pivot$Module <- factor(df_markersummary_pivot$Module, levels = colnames(scaled_expr), ordered = TRUE)
midpoint <- mean(df_markersummary_pivot$Expression)
marker_scalehm <- ggplot(df_markersummary_pivot, aes(x = Module, y = CellType, fill = Expression)) + geom_tile() + 
  theme_minimal() + ylab("") + xlab("") +
  scale_fill_gradient2(name = "Average\nScaled\nNormalized\nExpression", midpoint = midpoint,
                       low = "blue", mid = "white", high = "red") +
  scale_y_discrete(limits = rev) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.ticks = element_blank())
pdf("Marker_Modules_HM_AllCells_wCellType.pdf", width = 5, height = 4)
marker_scalehm
dev.off()

celltype_colorbar_df <- data.frame(CellType = levels(df_markersummary_pivot$CellType), X = 1)
celltype_colorbar_df$CellType <- factor(celltype_colorbar_df$CellType,
                                        levels = levels(df_markersummary_pivot$CellType), ordered = TRUE)
celltype_colorbar_ggbar <- ggplot(celltype_colorbar_df, aes(x = X, y = CellType, fill = CellType)) + geom_bar(position="fill", stat="identity") +
  theme_minimal() + scale_fill_manual(values = celltype_colors) + xlab("") + ylab("") + 
  scale_y_discrete(limits = rev) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_blank(), legend.position = "none")
pdf("Marker_Modules_HM_AllCells_wCellType.pdf", width = 6, height = 4)
plot_grid(celltype_colorbar_ggbar, marker_scalehm, align = "h", nrow = 1, rel_widths = c(1,3))
dev.off()

## Compare relative proportions of Histology (High- vs. Low-Grade) and Smoking Status (Current vs. Smoker) in each cell type
## Because of different numbers of nonimmune vs. immune plates sequenced per sample, comparison is done in nonimmune and immune cells separately
sce$Histology <- factor(sce$Histology, levels = unique(sce$Histology), ordered = FALSE)
sce$HistologyBinary <- ifelse(sce$Histology %in% c("Severe Dysplasia", "CIS", "LUSC"), "High Grade", "Low Grade")
sce$HistologyBinary <- factor(sce$HistologyBinary, levels = unique(sce$HistologyBinary), ordered = FALSE)
sce$Smoking_Status <- factor(sce$Smoking_Status, levels = unique(sce$Smoking_Status), ordered = FALSE)
sce$PCGA02_site <- factor(sce$PCGA02_site, levels = unique(sce$PCGA02_site), ordered = FALSE)
sce$Subject <- paste0("P", sce$Subject)
sce$Sample <- paste0("S", sce$Sample)
colData(sce) <- within(colData(sce), Histology <- relevel(Histology, ref = "Normal"))
colData(sce) <- within(colData(sce), HistologyBinary <- relevel(HistologyBinary, ref = "Low Grade"))
colData(sce) <- within(colData(sce), Smoking_Status <- relevel(Smoking_Status, ref = "Former"))
colData(sce) <- within(colData(sce), PCGA02_site <- relevel(PCGA02_site, ref = "Roswell"))
sce_nonimmune <- sce[,which(celdaClusters(sce) %in% c(1:18))]
sce_immune <- sce[,which(celdaClusters(sce) %in% c(19:30))]

cellshifts_binaryhistology_site_sccomp <- function(sce_object) {
  cat("Calculating Model with Site.\n")
  model_with_factor_association = 
    sce_object |>
    sccomp_estimate( 
      formula_composition = ~ 1 + HistologyBinary + Smoking_Status + PCGA02_site, 
      .sample = Sample, 
      .cell_group = CellType, 
      inference_method = "hmc",
      enable_loo = TRUE
    ) |> 
    sccomp_remove_outliers(cores = 4) |> 
    sccomp_test()
  
  plot_sum_histology_smoking <- model_with_factor_association |> 
    sccomp_boxplot(factor = "HistologyBinary")
  
  list(sitemodel = model_with_factor_association, summary_plots = plot_sum_histology_smoking)
}
nonimmune_binaryhistology_cellshifts_site <- cellshifts_binaryhistology_site_sccomp(sce_object = sce_nonimmune)
immune_binaryhistology_cellshifts_site <- cellshifts_binaryhistology_site_sccomp(sce_immune)

binary_histology_wsite <- rbind(nonimmune_binaryhistology_cellshifts_site$sitemodel, immune_binaryhistology_cellshifts_site$sitemodel)
binary_histology_wsite <- binary_histology_wsite %>% select(CellType:c_FDR)

## Figure 1D: Full model results bubble plot with stacked barplots of histology and smoking status distribution
df_hist <- data.frame(table(sce$CellType, sce$Histology))
df_hist$Var1 <- factor(df_hist$Var1,
                       levels = c("Basal Cells", "Goblet Cells", "Club Cells",
                                  "SMG Cells", "Ciliated Cells", "Fibroblasts", "Endothelial Cells", 
                                  "CD8+ T Cells", "CD4+ T Cells", "B Cells", "Plasma Cells", "Macrophages", 
                                  "Dendritic Cells", "Neutrophils", "Mast Cells"),
                       ordered = TRUE)
df_hist$Var2 <- factor(df_hist$Var2, levels = rev(levels(sce$Histology)), ordered = TRUE)
histology_sbp <- ggplot(df_hist, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") + coord_flip() + theme_classic() + 
  xlab("") + ylab("% Histology") + 
  scale_y_continuous(labels = scales::percent, expand = expansion(add = c(0,0.02))) + scale_x_discrete(expand = c(0,1), limits = rev) + scale_fill_manual(name = "Histology", values = hist_cols) +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank()) 

pdf("AllCells_CollapsedHistology_Boxplot.pdf", width = 5, height = 4)
histology_sbp
dev.off()

df_smoke <- data.frame(table(sce$CellType, sce$Smoking_Status))
df_smoke$Var1 <- factor(df_smoke$Var1,
                        levels = c("Basal Cells", "Goblet Cells", "Club Cells",
                                   "SMG Cells", "Ciliated Cells", "Fibroblasts", "Endothelial Cells", 
                                   "CD8+ T Cells", "CD4+ T Cells", "B Cells", "Plasma Cells", "Macrophages", 
                                   "Dendritic Cells", "Neutrophils", "Mast Cells"),
                        ordered = TRUE)
df_smoke$Var2 <- factor(df_smoke$Var2, levels = rev(levels(sce$Smoking_Status)), ordered = TRUE)
smoking_sbp <- ggplot(df_smoke, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") + coord_flip() + theme_classic() + 
  xlab("") + ylab("% Smoking Status") + 
  scale_y_continuous(labels = scales::percent, expand = expansion(add = c(0,0.02))) + scale_x_discrete(expand = c(0,1), limits = rev) + scale_fill_manual(name = "Smoking\nStatus", values = smoke_cols) +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank()) 
pdf("AllCells_Smoking_Boxplot.pdf", width = 5, height = 4)
smoking_sbp
dev.off()

stats_values <- binary_histology_wsite # For plotting, clip betas at [-3,3] and set low p-values to lowest value
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

stats_values <- stats_values[-which(stats_values$parameter %in% c("(Intercept)","PCGA02_siteUCL")),]
stats_values$parameter <- gsub("HistologyBinaryHigh Grade", "Severe Dysplasia - LUSC vs. Normal - Moderate Dysplasia", stats_values$parameter)
stats_values$parameter <- gsub("Smoking_StatusCurrent", "Current vs. Former Smoker", stats_values$parameter)
stats_values$parameter <- factor(stats_values$parameter,
                                 levels = c("Severe Dysplasia - LUSC vs. Normal - Moderate Dysplasia", "Current vs. Former Smoker"),
                                 ordered = TRUE)
effect_order <- order(unique(stats_values$c_effect[which(stats_values$parameter == "Severe Dysplasia - LUSC vs. Normal - Moderate Dysplasia")]))
stats_values$CellType <- factor(stats_values$CellType,
                                levels = rev(unique(stats_values$CellType)[effect_order]),
                                ordered = TRUE)

colors <- colorRampPalette(c("#053061", "#053061", "#053061", "#2166AC", "#4393C3", 
                                      "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#67001F", 
                                      "#67001F", "#67001F"))(200)
stats_values$CellTypeClusterOrder <- as.character(stats_values$CellType)
stats_values$CellTypeClusterOrder <- factor(stats_values$CellTypeClusterOrder,
                                            levels = c("Basal Cells", "Goblet Cells", "Club Cells",
                                                       "SMG Cells", "Ciliated Cells", "Fibroblasts", "Endothelial Cells", 
                                                       "CD8+ T Cells", "CD4+ T Cells", "B Cells", "Plasma Cells", "Macrophages", 
                                                       "Dendritic Cells", "Neutrophils", "Mast Cells"),
                                            ordered = TRUE)
cellshift_bubble_inclusterorder_noyaxis <- ggplot(stats_values, aes(x = parameter, y = CellTypeClusterOrder, label = sig, color = c_effect, 
                                             size = ifelse(is.na(logp), 2, logp), shape = ifelse(is.na(c_effect), "Missing", "Present"))) + 
  geom_point(pch=1, colour = "black", stroke = 2) +
  geom_point(stat = "identity", stroke = 2, na.rm = T, shape = 19) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 8),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 5),
        axis.line.y = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "", y = "", color = "Effect Size", size = "-log10(FDR)") + 
  guides(size = "legend", color = "colorbar") +
  geom_text(color = "white", size = 3, vjust = 0.9) +
  scale_x_discrete(labels = c("Severe Dysplasia - LUSC\nvs.\nNormal - Moderate Dysplasia", "Current vs.\nFormer Smoker")) +
  scale_y_discrete(limits = rev) +
  scale_shape_manual(values = c(Missing = 4, Present = 19)) +
  scale_color_gradientn(colours = colors, limits = c(-3, 3), breaks = c(-3, 0, 3)) 
pdf("SCComp_AllBiopsyCell_BubblePlot_BinaryHistology_040725.pdf", width = 5, height = 4)
cellshift_bubble_inclusterorder_noyaxis
dev.off()



########## ------------- SUPPLEMENTARY FIGURES/TABLES, FIGURE 1 ------------- ##########

## Supplementary Figure 1: Subject, Sample, TimePoint, Sample Type Dot Plot
sample_metadata <- openxlsx::read.xlsx("../Data_Management/PCGA2.0_scRNAseq_Biopsy+Brush_Metadata_2025-02-06_EEK.xlsx")
sample_metadata$Histology[which(sample_metadata$Histology %in% c("Denuded bronchial mucosa", "Inflammation", "Inflamation", "NAD"))] <- "Normal"
sample_metadata$Histology[which(sample_metadata$Histology %in% c("BASAL CELL HYPERPLASIA ", "Focal basal cell hyperplasia"))] <- "Basal Cell Hyperplasia"
sample_metadata$Histology[which(sample_metadata$Histology %in% c("METAPLASIA", "SqM"))] <- "Metaplasia"
sample_metadata$Histology[which(sample_metadata$Histology %in% c("MiD", "Mild dysplasia"))] <- "Mild Dysplasia"
sample_metadata$Histology[which(sample_metadata$Histology %in% c("MoD"))] <- "Moderate Dysplasia"
sample_metadata$Histology[which(sample_metadata$Histology %in% c("CIS+"))] <- "CIS"
sample_metadata$Histology[which(sample_metadata$Histology %in% c("INV"))] <- "LUSC"
sample_metadata$SampleType[which(sample_metadata$SampleType == "BronchialBrush")] <- "Bronchial Brush"
sample_metadata$SampleType[which(sample_metadata$SampleType == "NasalBrush")] <- "Nasal Brush"
sample_metadata$Anatomic_Site[which(is.na(sample_metadata$Anatomic_Site))] <- "Unknown"
sample_metadata$Subject[which(sample_metadata$Subject == "160 (EAR 005)")] <- "160"
sample_metadata$Subject[which(sample_metadata$Subject == "163 (EAR 009)")] <- "163"
sample_metadata$Subject[which(sample_metadata$Subject == "164 (EAR 010)")] <- "164"
sample_metadata$Anatomic_Site <- gsub(" ","", sample_metadata$Anatomic_Site)
sample_metadata$TimePoint <- as.character(sample_metadata$TimePoint)
sample_metadata$TimePoint <- factor(sample_metadata$TimePoint, 
                                    levels = as.character(1:8),
                                    ordered = TRUE)
num_plates <- data.frame(table(sample_metadata$PCGA02_ParentID))
colnames(num_plates) <- c("PCGA02_ParentID", "Num_Plates")
sample_tp_df <- sample_metadata %>% select(Subject, PCGA02_ParentID, TimePoint, Anatomic_Site, Histology, SampleType, PCGA02_site)
sample_tp_df <- full_join(x = sample_tp_df, y = num_plates)
sample_tp_df$Subject_Location = paste(sample_tp_df$Subject, sample_tp_df$Anatomic_Site)
sample_tp_df <- unique(sample_tp_df)

biopsy_tp_df <- sample_tp_df %>% filter(SampleType == "Biopsy")

# Match the rows of biopsies with brushes
for(i in seq(nrow(sample_tp_df))) {
  if(sample_tp_df[i, "SampleType"] %in% c("Bronchial Brush", "Nasal Brush")) {
    brush_subj <- sample_tp_df[i, "Subject"]
    biopsy_location <- which(sample_tp_df$Subject == brush_subj & sample_tp_df$SampleType == "Biopsy")
    if(length(biopsy_location) > 0) { # Some subjects don't have associated biopsies
      m <- biopsy_location[1]
      sample_tp_df[i, "Subject_Location"] <- sample_tp_df[m, "Subject_Location"]
    } else { 
      sample_tp_df[i, "Subject_Location"] <- sample_tp_df[i, "Subject"]
    }
  }
}

sample_tp_df$Subject <- as.numeric(sample_tp_df$Subject)
sample_tp_df <- sample_tp_df %>% arrange(desc(Subject))
sample_tp_df$Subject <- as.character(sample_tp_df$Subject)
sample_tp_df$Subject_Location <- factor(sample_tp_df$Subject_Location,
                                        levels = unique(sample_tp_df$Subject_Location),
                                        ordered = TRUE) # orders x-axis smallest to largest by numeric subject number
sample_tp_df$Histology <- factor(sample_tp_df$Histology,
                                 levels = c("Normal", "Basal Cell Hyperplasia", "Metaplasia", "Mild Dysplasia",
                                            "Moderate Dysplasia", "Severe Dysplasia", "CIS", "LUSC", 
                                            "Brush NO HISTOLOGY ", "dysplastic squamous epithelium "),
                                 ordered = TRUE)

histology_colors <- readRDS("histology_colors.rds")
histology_colors[["Brush NO HISTOLOGY "]] <- "black"
sample_timepoint_dotplot <- ggplot(sample_tp_df, aes(x = TimePoint, y = Subject_Location, color = Histology, size = Num_Plates)) + geom_point(shape = 15) +
  theme_classic() + facet_wrap(~SampleType) + scale_size_continuous(name = "Number of Plates") +
  scale_color_manual(values = histology_colors) + xlab("Time Point") + ylab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
sample_timepoint_dotplot

sample_tp_df$X <- 1
site_cols <- c("UCL" = "purple1", "Roswell" = "gold")
site_label_gg <- ggplot(sample_tp_df, aes(x = X, y = Subject_Location, fill = PCGA02_site)) + geom_bar(position="fill", stat="identity") +
  theme_minimal() + xlab("") + ylab("Subject, Anatomic Location") + scale_fill_manual(name = "Collection Site", values = site_cols) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_blank())
site_label_gg

pdf("Subject_TimePoint_Histology_SampleType_SummaryPlot.pdf", width = 10, height = 10)
plot_grid(plotlist = list(site_label_gg, sample_timepoint_dotplot), axis = "tblr", align = "hv", nrow = 1, rel_widths = c(1,2))
dev.off()

## Supplementary Figure 2: Cluster UMAP
umap_coords <- reducedDim(altExp(sce), "celda_UMAP")
colnames(umap_coords) <- c("UMAP1", "UMAP2")
cluster_df <- cbind(data.frame(Clusters = paste0("K", celdaClusters(sce))), umap_coords)
cluster_df$Clusters <- factor(cluster_df$Clusters, levels = paste0("K",levels(celdaClusters(sce))), ordered = TRUE)
cluster_celltype_df <- cbind(data.frame(Clusters = paste0(sce$CellType, " (K", celdaClusters(sce),")")), umap_coords)
celltype_df <- cbind(data.frame(CellType = sce$CellType), umap_coords)
cluster_colors <- distinctColors(
  n = length(unique(celdaClusters(sce))),
  hues = c("red", "cyan", "orange", "pink1", "blue", "yellow", "purple", "green", "forestgreen", "magenta"))
cluster_colors2 <- cluster_colors
cluster_colors2[6] <- cluster_colors[9]; cluster_colors2[9] <- cluster_colors[6]
cluster_colors2[24] <- cluster_colors[29]; cluster_colors2[29] <- cluster_colors[24]
cluster_colors2[1] <- cluster_colors[28]; cluster_colors2[28] <- cluster_colors[1]
saveRDS(cluster_colors2, "/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/all_cells_cluster_colors.rds")

cluster_umap <- ggplot(cluster_df, aes(x = UMAP1, y = UMAP2, color = Clusters)) + geom_point(size = 0.25) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                          axis.ticks = element_blank()) + scale_color_manual(values = cluster_colors2)
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
                          axis.ticks = element_blank()) + scale_color_manual(values = cluster_colors2)
cluster_celltype_umap_labeled <- add_labels(cluster_celltype_df, "Clusters", plt = cluster_celltype_umap)
cluster_celltype_umap_labeled <- cluster_celltype_umap_labeled + labs(color = "Cluster")
pdf("AllCell_CellType+Cluster_UMAP.pdf", width = 8, height = 6)
cluster_celltype_umap_labeled
dev.off()

## Supplementary Figure 3: Module Probability Heatmap
mod_map <- celdaProbabilityMap(sce, useAssay = "decontXcounts")
pdf("AllBiopsyCell_ProbabilityMap.pdf", width = 10, height = 10)
mod_map
dev.off()

## Supplementary Table 4: Enrichment Terms for All Biopsy Cell Modules (EnrichR)
modules <- read.table("070124_celda_K30_L104.tsv", sep = "\t", header = TRUE)
modules <- lapply(as.list(modules), function(x) {
  x <- x[x != ""]
  x
})
modules_enrichr <- lapply(modules, function(x) {
  enrichments <- enrichr(genes = x, databases = "GO_Biological_Process_2023")
  enrichments
})
modules_enrichr_df <- data.frame(Module = paste0("L",1:104), 'Top 10 Enriched Genes' = rep("",104), 'Top 5 Enrichment Terms' = rep("",104), check.names = FALSE)
for(i in seq_along(modules_enrichr_df$Module)) {
  df <- modules_enrichr[[i]][["GO_Biological_Process_2023"]]
  genes <- paste(modules_enrichr[[i]][["GO_Biological_Process_2023"]][["Genes"]][1:5], collapse = ";")
  genes <- strsplit(genes, split = ";")[[1]]
  genes <- unique(genes)
  genes <- genes[1:min(length(genes),10)]
  modules_enrichr_df$`Top 10 Enriched Genes`[i] <- paste(genes, collapse = ";")
  modules_enrichr_df$`Top 5 Enrichment Terms`[i] <- paste(modules_enrichr[[i]][["GO_Biological_Process_2023"]][["Term"]][1:5], collapse = ";")
}
write.xlsx(modules_enrichr_df, file = "Biopsy_Module_Enrichr.xlsx")

## Supplementary Table 5: Cell Type x Histology and Smoking Status (SCComp)
write.xlsx(binary_histology_wsite, file = "CellType_Immune+Nonimmune_wSite_FullModel_040725.xlsx", rowNames = FALSE) 

saveRDS(sce, "AllBiopsyCell_SCE_Final.rds")
