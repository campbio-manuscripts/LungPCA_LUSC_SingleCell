library(celda)
library(singleCellTK)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)
library(scater)
library(nlme)
library(glmmTMB)
library(ggplot2)
library(sccomp)
library(scales)
library(gridExtra)
library(cowplot)
library(patchwork)
library(stringr)
library(tidyr)

##############################################################################################################
##################################### Read in data

# sce$immune_CellType[which(sce$immune_celda_cluster == 2)] <- "Eosinophils"
# altExp(sce)$immune_CellType <- sce$immune_CellType
# saveRDS(sce, "./040725_PCGA_ImmuneCells_NoMTPseudo_K27_L57_annotated_altExpIGgenes.rds")

sce <- readRDS("./040725_PCGA_ImmuneCells_NoMTPseudo_K27_L57_annotated_altExpIGgenes.rds")
sce_orig <- sce

# celltype_colors <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/celltype_colors.rds")
# celltype_colors <- c(celltype_colors, "#0000C1", "#009EAF")
# names(celltype_colors)[33] <- "Proliferating T Cells"
# names(celltype_colors)[34] <- "Eosinophils"
# saveRDS(celltype_colors, "../celltype_colors.rds")

celltype_colors <- readRDS("../celltype_colors.rds")

##############################################################################################################
##################################### Figure 4A
reducedDim(sce, "celda_UMAP") <- reducedDim(altExp(sce), "celda_UMAP")

pdf("./Figures/Fig4A_Immune_CellType_UMAP.pdf", width = 13, height = 9)
plotSCEDimReduceColData(sce, 
                        colorBy = "immune_CellType", 
                        reducedDimName = "celda_UMAP", 
                        colorScale = celltype_colors,
                        dotSize = .7, 
                        labelClusters = FALSE) + 
  theme_classic()
dev.off()

##############################################################################################################
##################################### Figure 4B
################### Organize object
f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")), 
                                    colData = colData(sce))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")

cellprobs <- data.frame(t(assay(sce_modular,"LogNormalize")))
cellprobs$Cluster <- sce_modular$immune_CellType ## 

mean_cellprobs <- cellprobs %>% group_by(Cluster) %>% summarize_all(mean)
mean_cellprobs <- data.frame(t(mean_cellprobs))
colnames(mean_cellprobs) <- mean_cellprobs["Cluster",]
mean_cellprobs <- mean_cellprobs[2:56,] ## 
mean_cellprobs <- sapply(mean_cellprobs, as.numeric)
mean_cellprobs <- data.frame(mean_cellprobs, row.names = colnames(cellprobs)[1:55])
cellprobs$Sample <- sce$Sample
cellprobs$Smoking_Status <- sce$Smoking_Status

################### Find modules
immune_clust_de <- function(expr = cellprobs, group_var = "Cluster", clusters) {
  # expr_subset <- expr[which(expr[,group_var] %in% clusters),]
  expr_subset <- expr
  each_model_res <- lapply(clusters, function(x) {
    expr_subset[,"ImmuneCluster"] <- ifelse(expr_subset[,group_var] == x, x, paste0("Not",x))
    expr_subset[,"ImmuneCluster"] <- factor(expr_subset[,"ImmuneCluster"],
                                            levels = c(paste0("Not",x),x),
                                            ordered = TRUE)
    model_p <- data.frame(Cluster = x, 
                          Module = names(expr_subset)[1:57], 
                          P = rep(0, 57),
                          Beta = rep(0, 57), 
                          Group_Mean = rep(0, 57), 
                          NotGroup_Mean = rep(0, 57), 
                          row.names = names(expr_subset)[1:57]) # Smoking_P = rep(0, 55), Smoking_Beta = rep(0, 55), 
    for(i in rownames(model_p)) {
      model_p[i,"Group_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] == x),i], na.rm = TRUE)
      model_p[i,"NotGroup_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] != x),i], na.rm = TRUE)
      
      model_construct <- as.formula(paste0(i," ~ ImmuneCluster"))
      # if(group_var == "Smoking_Status") {
      #   model_construct <- as.formula(paste0(i," ~ ImmuneCluster"))
      # } else {
      #   model_construct<- as.formula(paste0(i," ~ ImmuneCluster + Smoking_Status"))
      # }
      model <- lme(model_construct,
                   random = ~ 1|Sample,
                   data = expr_subset,
                   na.action = na.omit)
      tTab <- summary(model)$tTable
      model_p[i,"Beta"] <- tTab["ImmuneCluster.L","Value"]
      model_p[i,"P"] <- tTab["ImmuneCluster.L","p-value"]
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

diff_modules <- immune_clust_de(clusters = unique(cellprobs$Cluster))
diff_modules_res <- Reduce(diff_modules, f = rbind)
diff_modules_res_significant <- diff_modules_res[which(diff_modules_res$SignificantPositive),]

saveRDS(diff_modules_res_significant, "./Objects_Results/Fig4B_diffExpModules_byCellType_Significant.rds")

################### Specify order of cell types
celltype_order <- c("CD4+ T Cells", "CD8+ T Cells", "Proliferating T Cells", "NK Cells", "B Cells", "Dendritic Cells", "Plasma Cells",  
                    "Plasmacytoid Dendritic Cells", "Monocytes", "Macrophages", "Neutrophils", "Eosinophils", "Mast Cells", "L27+ Cells")
celltype_modules <- as.data.frame(matrix(NA, ncol = 2, nrow = length(unique(sce$immune_CellType))))
celltype_modules$V1 <- celltype_order

for (i in 1:nrow(celltype_modules)) {
  celltype <- diff_modules_res_significant[which(diff_modules_res_significant$Cluster == celltype_modules$V1[i]),]
  celltype_modules$V2[i] <- celltype$Module[which(celltype$Beta == max(celltype$Beta))]
}

################### Prep heatmap data
m <- unique(celltype_modules$V2)
m <- gsub("L23", "L21", m)

scaled_expr <- t(assay(sce_modular,"LogNormalize")[m,])
scaled_expr <- scale(scaled_expr)

colnames(scaled_expr)[which(colnames(scaled_expr) == "L46")] <- "L46 (IL32, CD3D)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L45")] <- "L45 (CD8A)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L49")] <- "L49 (MKI67)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L29")] <- "L29 (GNLY, NKG7)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L54")] <- "L54 (MS4A1)"
# colnames(scaled_expr)[which(colnames(scaled_expr) == "L23")] <- "L23 (MS4A6A)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L21")] <- "L21 (FCER1A)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L53")] <- "L53 (MZB1)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L18")] <- "L18 (IRF8)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L22")] <- "L22 (LYZ)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L19")] <- "L19 (C1QA)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L15")] <- "L15 (CSF3R)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L56")] <- "L56 (CLC)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L55")] <- "L55 (CPA3, MS4A2)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L27")] <- "L27 (LRRK2, APOE)"

df <- cbind(data.frame(CellType = sce_modular$immune_celda_cluster), scaled_expr)
df_markersummary <- df %>% group_by(CellType) %>% summarize_all(mean)
df_markersummary_pivot <- df_markersummary %>% pivot_longer(!CellType, names_to = "Module", values_to = "Expression")
df_markersummary_pivot$CellType <- paste0("K", df_markersummary_pivot$CellType)
df_markersummary_pivot$CellType <- factor(df_markersummary_pivot$CellType,
                                          levels = c("K14", "K5", "K2", "K3", "K1", "K4",
                                                     "K27", "K26", "K24", "K25", "K15", 
                                                     "K21", "K20", "K22", "K23", "K7", "K13",
                                                     "K17", "K18", "K8", "K12", "K6", "K16",
                                                     "K11", "K19", "K10", "K9"),
                                          ordered = TRUE)

df_markersummary_pivot$Module <- factor(df_markersummary_pivot$Module, levels = colnames(scaled_expr), ordered = TRUE)
midpoint <- mean(df_markersummary_pivot$Expression)

df_markersummary_pivot$Expression[which(df_markersummary_pivot$Expression > 2)] <- 2

################### Plot heatmap
marker_scalehm <- ggplot(df_markersummary_pivot, aes(x = Module, y = CellType, fill = Expression)) + geom_tile() + 
  theme_minimal() + ylab("") +
  scale_fill_gradient2(name = "Average\nScaled\nNormalized\nExpression", midpoint = midpoint,
                       low = "blue", mid = "white", high = "red") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        # axis.text.y = element_blank(), 
        # plot.margin = unit(c(0, 0, 0, -0.5), "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.ticks = element_blank())

pdf("./Figures/Fig4B_CellType_Modules_Heatmap_temp.pdf", width = 4, height = 7)
marker_scalehm
dev.off()

##############################################################################################################
##################################### Figure 4C

sce <- sce_orig

################### Remove any cells that don't have a low/mixed/high grade label
sce <- sce[,which(sce$BasalGroupMixed %in% c("Low Grade Basal Sample", "High Grade Basal Sample"))]
sce$BasalGroupMixed <- gsub("Low Grade Basal Sample", "LowGradeBasalSample", sce$BasalGroupMixed)
sce$BasalGroupMixed <- gsub("High Grade Basal Sample", "HighGradeBasalSample", sce$BasalGroupMixed)
sce$BasalGroupMixed <- factor(sce$BasalGroupMixed, levels = c("LowGradeBasalSample", "HighGradeBasalSample"))

sce$Smoking_Status <- unfactor(sce$Smoking_Status)
sce$Smoking_Status <- factor(sce$Smoking_Status, levels = c("Former", "Current"))
sce$PCGA02_site <- factor(sce$PCGA02_site, levels = c("Roswell", "UCL"))
sce$immune_celda_cluster <- paste0("K", sce$immune_celda_cluster, sep = "")
sce$Sample <- paste0("S", sce$Sample, sep = "")

################### Run SCComp on celda cluster
model_with_factor_association = 
  sce |>
  sccomp_estimate(  #sccomp_glm
    formula_composition = ~ 1 + BasalGroupMixed + Smoking_Status + PCGA02_site, 
    # formula_variability = ~ 1,
    .sample = Sample, 
    .cell_group = immune_celda_cluster, 
    # bimodal_mean_variability_association = TRUE,
    # cores = 1, 
    enable_loo = TRUE,
    inference_method = "hmc"
    # variational_inference = FALSE
  ) |> 
  sccomp_remove_outliers(cores = 4) |>
  sccomp_test()

saveRDS(model_with_factor_association, "./Objects_Results/sccomp_results_BasalGroupMixed_covSmokingSite.rds")

################### Organize the summary table
stats_values <- model_with_factor_association

stats_values$c_FDR <- as.numeric(as.vector(stats_values$c_FDR))
stats_values$c_FDR[stats_values$c_FDR == 0] <- 2e-16
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

stats_values$immune_celda_cluster <- factor(stats_values$immune_celda_cluster, levels = c("K14", "K5", "K2", "K3", "K1", "K4",
                                                                                          "K27", "K26", "K24", "K25", "K15", 
                                                                                          "K21", "K20", "K22", "K23", "K7", "K13",
                                                                                          "K17", "K18", "K8", "K12", "K6", "K16",
                                                                                          "K19", "K10", "K9", "K11"))

stats_values <- stats_values[-which(stats_values$parameter == "(Intercept)"),]
stats_values <- stats_values[-which(stats_values$parameter == "PCGA02_siteUCL"),]
stats_values$parameter <- gsub("BasalGroupMixedHighGradeBasalSample", "HGB vs LGB", stats_values$parameter)
stats_values$parameter <- gsub("Smoking_StatusCurrent", "SmCurr", stats_values$parameter)
stats_values$parameter <- factor(stats_values$parameter, levels = c("HGB vs LGB", "SmCurr"))

################### Plot
########## Bubbleplot
colors <- colorRampPalette(c("#053061", "#053061", "#053061", "#2166AC", "#4393C3", 
                             "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#67001F", 
                             "#67001F", "#67001F"))(200)

rescales <- rescale(c(min(na.omit(stats_values$c_effect)), 0, max(na.omit(stats_values$c_effect))), 
                    limits = c(min(na.omit(stats_values$c_effect)), max(na.omit(stats_values$c_effect))))
p <- ggplot(stats_values, aes(x = parameter, y = immune_celda_cluster, label = sig, color = c_effect, 
                              size = ifelse(is.na(logp), 2, logp), shape = ifelse(is.na(c_effect), "Missing", "Present"))) + 
  geom_point(pch=1, colour = "black", stroke = 2) +
  geom_point(stat = "identity", stroke = 2, na.rm = T, shape = 19) + 
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5)) +
  guides(color="none") +
  labs(x = "", y = "", color = "Effect Size") + 
  guides(size = F, shape = F) +
  geom_text(color = "white", size = 3, vjust = 0.9) +
  scale_shape_manual(values = c(Missing = 4, Present = 19)) +
  scale_color_gradientn(colours = colors, values = rescales)

########## Boxplots
cd <- as.data.frame(colData(sce))
cd$immune_celda_cluster <- factor(cd$immune_celda_cluster, levels = c("K14", "K5", "K2", "K3", "K1", "K4",
                                                                      "K27", "K26", "K24", "K25", "K15", 
                                                                      "K21", "K20", "K22", "K23", "K7", "K13",
                                                                      "K17", "K18", "K8", "K12", "K6", "K16",
                                                                      "K19", "K10", "K9", "K11"))

cd$BasalGroup_rename <- cd$BasalGroupMixed
cd$BasalGroup_rename <- gsub("HighGradeBasalSample", "HG Basal Sample", cd$BasalGroup_rename)
cd$BasalGroup_rename <- gsub("LowGradeBasalSample", "LG Basal Sample", cd$BasalGroup_rename)
cd$BasalGroup_rename <- ordered(cd$BasalGroup_rename, levels = c("HG Basal Sample", "LG Basal Sample"))

cd_basal <- cd %>% dplyr::group_by(immune_celda_cluster) %>% dplyr::count(BasalGroup_rename) %>% na.omit()

basal <- ggplot(cd_basal, aes(y = immune_celda_cluster, fill = BasalGroup_rename, x = n)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  scale_fill_manual(values = c("LG Basal Sample" = "blue",
                               # "Mixed Grade Basal Sample" = "#E0B0FF",
                               "HG Basal Sample" = "red"),
                    name = "Type of Basal Cells in Sample") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(fill = "none") +
  labs(x = "")  +
  scale_x_continuous(labels = scales::percent)

cd_smoke <- cd %>% dplyr::group_by(immune_celda_cluster) %>% dplyr::count(Smoking_Status) %>% na.omit()

smoke <- ggplot(cd_smoke, aes(y = immune_celda_cluster, fill = Smoking_Status, x = n)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  scale_fill_manual(values = c("Former" = "grey40", "Current" = "red")) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(fill = "none") +
  labs(x = "") +
  scale_x_continuous(labels = scales::percent)

################### Combine & save plots
pdf(paste0("./Figures/Fig4C_bubbleplot_BasalGroupMixed_covSmokingSite.pdf"), height = 8, width = 6)
grid.arrange(p, smoke, basal, ncol = 3, nrow=1, widths = c(1, 0.5, 0.5))
dev.off()


##############################################################################################################
##################################### Figure 4D-F

########################## Figure 4D-F: Boxplots
percent_celltype <- data.frame(prop.table(table(sce$Sample, sce$immune_celda_cluster), margin = 1))
colnames(percent_celltype) <- c("Sample", "KCluster", "Proportion")
sample_metadata_df <- data.frame(unique(data.frame("Sample" = sce$Sample, "SmokingStatus" = sce$Smoking_Status, "BasalGroup" = sce$BasalGroupMixed)))
percent_celltype <- left_join(percent_celltype, 
                              sample_metadata_df,
                              by = c("Sample"))
cluster_metadata_df <- data.frame(unique(data.frame("KCluster" = sce$immune_celda_cluster, "CellType" = sce$immune_CellType)))
percent_celltype <- left_join(percent_celltype, 
                              cluster_metadata_df,
                              by = c("KCluster"))
percent_celltype$KCluster_CellType <- paste(percent_celltype$KCluster, percent_celltype$CellType, sep = " - ")
percent_celltype$BasalGroup <- factor(percent_celltype$BasalGroup,
                                      levels = c("Low Grade Basal Sample", "High Grade Basal Sample"),
                                      ordered = TRUE)
################### Plot boxplots
pdf("./Figures/Fig4D_K26_Sccomp_Boxplot.pdf", height = 2)
df <- percent_celltype[which(percent_celltype$KCluster == "26"),]
df <- df[-which(is.na(df$BasalGroup)),]
ggplot(df, aes(x = BasalGroup, y = Proportion, fill = BasalGroup)) + geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(labels = c("Low Grade\nBasal Sample", "High Grade\nBasal Sample")) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Low Grade Basal Sample" = "blue", "High Grade Basal Sample" = "red")) +
  labs(title = "K26", y = "% Cells", x = "") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  coord_flip()
dev.off()

pdf("./Figures/Fig4E_K3_Sccomp_Boxplot.pdf", height = 2)
df <- percent_celltype[which(percent_celltype$KCluster == "3"),]
df <- df[-which(is.na(df$BasalGroup)),]
ggplot(df, aes(x = BasalGroup, y = Proportion, fill = BasalGroup)) + geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(labels = c("Low Grade\nBasal Sample", "High Grade\nBasal Sample")) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Low Grade Basal Sample" = "blue", "High Grade Basal Sample" = "red")) +
  labs(title = "K3", y = "% Cells", x = "") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  coord_flip()
dev.off()

pdf("./Figures/Fig4E_K4_Sccomp_Boxplot.pdf", height = 2)
df <- percent_celltype[which(percent_celltype$KCluster == "4"),]
df <- df[-which(is.na(df$BasalGroup)),]
ggplot(df, aes(x = BasalGroup, y = Proportion, fill = BasalGroup)) + geom_boxplot() +
  theme_classic() +
  scale_x_discrete(labels = c("Low Grade\nBasal Sample", "High Grade\nBasal Sample")) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Low Grade Basal Sample" = "blue", "High Grade Basal Sample" = "red")) +
  labs(title = "K4", y = "% Cells", x = "") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  coord_flip()
dev.off()

pdf("./Figures/Fig4F_K13_Sccomp_Boxplot.pdf", height = 2)
df <- percent_celltype[which(percent_celltype$KCluster == "13"),]
df <- df[-which(is.na(df$BasalGroup)),]
ggplot(df, aes(x = BasalGroup, y = Proportion, fill = BasalGroup)) + geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(labels = c("Low Grade\nBasal Sample", "High Grade\nBasal Sample")) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Low Grade Basal Sample" = "blue", "High Grade Basal Sample" = "red")) +
  labs(title = "K13", y = "% Cells", x = "") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  coord_flip()
dev.off()

########################## Figure 4D-F: Heatmaps
sce <- sce_orig

########## Read in cell type signatures
wu_celltype_signatures <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/wu_celltype_signatures.rds")
xie_neutrophil_geneLists <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/xie_neutrophil_geneLists.rds")
names(xie_neutrophil_geneLists) <- paste0("Xie_", names(xie_neutrophil_geneLists), sep = "")
panglao_celltype_modules <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/panglao_celltype_modules.rds")
names(panglao_celltype_modules) <- paste0("Panglao_", names(panglao_celltype_modules), sep = "")
sikkema_cellType_markers <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/sikkema_cellType_markers.rds")
names(sikkema_cellType_markers) <- paste0("Sikkema_", names(sikkema_cellType_markers), sep = "")
hao_bCell_plasmaCell_signatures <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/hao_bCell_plasmaCell_signatures.rds")
names(hao_bCell_plasmaCell_signatures) <- paste0("Hao_", names(hao_bCell_plasmaCell_signatures), sep = "")
sinjab_celltype_modules <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/sinjab_celltype_modules.rds")
names(sinjab_celltype_modules) <- paste0("Sinjab_", names(sinjab_celltype_modules), sep = "")
travaglini_celltype_modules <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/travaglini_celltype_modules.rds")
Mascaux2019_Supp4_ImmuneSigs <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/040725_Mascaux2019_Supp4_ImmuneSigs.rds")
names(Mascaux2019_Supp4_ImmuneSigs) <- paste0("Mascaux_", names(Mascaux2019_Supp4_ImmuneSigs), sep = "")
science_TIB_bcell <- readRDS("/restricted/projectnb/pcga/TargetedBCR/DECAMP_snuc_TRUST4/ExternalData/2024Science_Pancancer_TIB/PanB_gene_sets.rds")
names(science_TIB_bcell) <- paste0("Science_BCell_", names(science_TIB_bcell), sep = "")
science_TIB_plasma <- readRDS("/restricted/projectnb/pcga/TargetedBCR/DECAMP_snuc_TRUST4/ExternalData/2024Science_Pancancer_TIB/PanB_PCsubset_DE_gene_sets.rds")
names(science_TIB_plasma) <- paste0("Science_Plasma_", names(science_TIB_plasma), sep = "")
cell_bcell_plasma <- readRDS("/restricted/projectnb/pcga/TargetedBCR/DECAMP_snuc_TRUST4/ExternalData/2024Cell_Pancancer_TIB/Bsubset_curated_gene_signatures.rds")
names(cell_bcell_plasma) <- paste0("Cell_", names(cell_bcell_plasma), sep = "")
spatial_markers <- readRDS("/restricted/projectnb/pcga/TargetedBCR/DECAMP_snuc_TRUST4/ExternalData/2023NG_Spatial_LungCellAtlas/Immune_cell_markers.rds")
names(spatial_markers) <- paste0("Madissoon_", names(spatial_markers), sep = "")
spatial_markers$Madissoon_Macro_CCL <- unique(spatial_markers$Madissoon_Macro_CCL)
zheng_modules <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/zheng_modules.rds")
names(zheng_modules) <- paste0("Zheng_", names(zheng_modules), sep = "")
dykema_list <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/dykema_tRegCell_signatures.rds")
names(dykema_list) <- paste0("Dykema_", names(dykema_list), sep = "")
zilionis_neut <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/Zilionis_Fig3_NeutrophilSigs.rds")
names(zilionis_neut) <- paste0("Zilionis_", names(zilionis_neut), sep = "")
zilionis_momac <- readRDS("/restricted/projectnb/pcga/RConrad/Imported_datasets/Zilionis_Fig5_MoMacSigs.rds")
names(zilionis_momac) <- paste0("Zilionis_", names(zilionis_momac), sep = "")

nan_5B <- c("RBP7", "PADI4", "S100A12", "MME", "APMAP", "TMX4", "ADAM19", "coro1A", "VNN2", "PROK2", "CDA", "MMP9", "FLOT2", "SELL",
            "MSRB1", "CXCR1", "UBN1", "S100A8", "S100A9", "TKT", "LRP10", "PLXNC1", "MEGF9", "LRRK2", "IL17RA", "CNN2", "S100A4", 
            "CMTM2", "RTN3", "GCA")
tan_5B <- c("RPL23", "SQSTM1", "CSTB", "EEF1A1", "PTMA", "RPS28", "NAP1L1", "DYNLL1", "PLEKHB2", "HSP90AB1", "CD74", "ST13", "CXCL2",
            "CSF1", "BHLHE40", "TNF", "COX20", "OLR1", "IRAK2", "CD22", "HLA-DRA", "CDKN1A", "CCRL2", "LGALS3", "CD83", "TGM2", "CCL4",
            "CCL4L2", "C15orf48", "CCL3")
nan_6B <- c("CXCR1", "CXCR2", "PTGS2", "SELL", "CSF3R", "FCGR3B")
tan_6B <- c("OLR1", "VEGFA", "CD83", "ICAM1", "CXCR4")
nan_s4 <- c("AGO4", "ARG1", "CYP4F3", "ERGIC1", "FRAT2", "LRP10", "MGAM", "MMP25", "MSRB1",
            "NDEL1", "NFE2", "PADI4", "PBX2", "PHOSPHO1", "RASGRP4", "REPS2", "SULT1B1", "TSEN34", "XKR8")
tan_s4 <- c("CCR3", "CCRL2", "DDIT3", "FLOT1", "HIF1A", "IRAK2", "MAFF", "MAP1LC3B2", "MCOLN1",
            "NBN", "NOD2", "PI3", "PLAU", "PPIF", "TGM3", "TOM1", "UBR5-AS1", "ZNF267")

nan_tan <- list(nan_5B, tan_5B, nan_6B, tan_6B, nan_s4, tan_s4)
names(nan_tan) <- c("NAN_5B", "TAN_5B", "NAN_6B", "TAN_6B", "NAN_S4", "TAN_S4")


################################################################################
########## Calculate VAM scores for all signatures
sigs <- c(wu_celltype_signatures, xie_neutrophil_geneLists, panglao_celltype_modules, 
          sikkema_cellType_markers, hao_bCell_plasmaCell_signatures, sinjab_celltype_modules, 
          travaglini_celltype_modules, Mascaux2019_Supp4_ImmuneSigs, science_TIB_bcell, science_TIB_plasma, cell_bcell_plasma, spatial_markers,
          zheng_modules, dykema_list, zilionis_neut, zilionis_momac, nan_tan) 

# sigs <- c(madissoon_immune)

names(sigs) <- gsub(" ", "_", names(sigs))

########## Implement expression filter & score gene sets
sce_scoring <- sce
sce_scoring <- selectFeatures(sce_scoring, useAssay = "decontXcounts") 
sigs_filtered <- lapply(sigs, function(x) {
  x <- x[x %in% rowData(altExp(sce_scoring))$Gene_name]
})
sigs_filtered2 <- Filter(function(x) length(x) > 0, sigs_filtered)

altExp(sce_scoring) <- importGeneSetsFromList(altExp(sce_scoring), geneSetList = sigs_filtered2, collectionName = "GeneSetCollection", by = "Gene_name")
altExp(sce_scoring) <- runVAM(inSCE = altExp(sce_scoring), geneSetCollectionName = "GeneSetCollection", useAssay = "LogNormalize")

scores <- as.data.frame(altExp(sce_scoring)@int_colData@listData[["reducedDims"]]@listData[["VAM_GeneSetCollection_CDF"]])

saveRDS(scores, "./Objects_Results/Fig4DEF_External_CellTypeSignatures_VAM_Scores.rds")

########## Plot VAM scores for all signatures across all clusters via violin plot
scores <- readRDS("./Objects_Results/Fig4DEF_External_CellTypeSignatures_VAM_Scores.rds")

colData(sce) <- cbind(colData(sce), scores)
cd <- as.data.frame(colData(sce))

colnames(cd) <- gsub("\\.", "", colnames(cd))
colnames(scores) <- gsub("\\.", "", colnames(scores))
colnames(cd) <- gsub("\\(", "", colnames(cd))
colnames(scores) <- gsub("\\(", "", colnames(scores))
colnames(cd) <- gsub("\\)", "", colnames(cd))
colnames(scores) <- gsub("\\)", "", colnames(scores))
colnames(cd) <- gsub("\\-", "", colnames(cd))
colnames(scores) <- gsub("\\-", "", colnames(scores))
colnames(cd) <- gsub("\\+", "", colnames(cd))
colnames(scores) <- gsub("\\+", "", colnames(scores))
colnames(cd) <- gsub("\\/", "", colnames(cd))
colnames(scores) <- gsub("\\/", "", colnames(scores))
colnames(cd) <- gsub("\\Ü", "U", colnames(cd))
colnames(scores) <- gsub("\\Ü", "U", colnames(scores))
colnames(cd) <- gsub("\\œ", "o", colnames(cd))
colnames(scores) <- gsub("\\œ", "o", colnames(scores))
colnames(cd) <- gsub("\\ï", "i", colnames(cd))
colnames(scores) <- gsub("\\ï", "i", colnames(scores))
colnames(cd) <- gsub("\\:", "", colnames(cd))
colnames(scores) <- gsub("\\:", "", colnames(scores))
colnames(cd) <- gsub("\\'", "", colnames(cd))
colnames(scores) <- gsub("\\'", "", colnames(scores))
colnames(cd) <- gsub("\\&", "", colnames(cd))
colnames(scores) <- gsub("\\&", "", colnames(scores))

pdf("./Figures/Fig4DEF_External_CellTypeSignatures_VAM_Scores_ViolinPlots.pdf", width = 15, height = 8)
for (i in colnames(scores)) {
  print(ggplot(data = cd, aes_string(x = "immune_celda_cluster", y = i, fill = "immune_CellType")) +
          geom_violin() + 
          theme_classic() +
          scale_fill_manual(values = celltype_colors) +
          ggtitle(i)
  )
}
dev.off()

################################################################################
########################## Identify modules differentially expressed between clusters of the same cell type
f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")), 
                                    colData = colData(sce))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")

cellprobs <- data.frame(t(assay(sce_modular,"LogNormalize")))
cellprobs$Cluster <- sce_modular$immune_celda_cluster ## 

cellprobs$Sample <- sce$Sample
cellprobs$Smoking_Status <- sce$Smoking_Status
cellprobs$Site <- sce$PCGA02_site

immune_clust_de <- function(expr = cellprobs, group_var = "Cluster", clusters, beta = 0.5) {
  # expr_subset <- expr[which(expr[,group_var] %in% clusters),]
  expr_subset <- expr
  each_model_res <- lapply(clusters, function(x) {
    # expr_subset[,"ImmuneCluster"] <- ifelse(expr_subset[,group_var] == x, x, paste0("Not",x))
    expr_subset[,"ImmuneCluster"] <- ifelse(expr_subset[,group_var] == x, "B", "A")
    # expr_subset[,"ImmuneCluster"] <- factor(expr_subset[,"ImmuneCluster"],
    #                                         levels = c(paste0("Not",x),x),
    #                                         ordered = TRUE)
    model_p <- data.frame(Cluster = x, 
                          Module = names(expr_subset)[1:57], 
                          P = rep(0, 57), 
                          Beta = rep(0, 57), 
                          Group_Mean = rep(0, 57), 
                          NotGroup_Mean = rep(0, 57), 
                          row.names = names(expr_subset)[1:57]) # Smoking_P = rep(0, 55), Smoking_Beta = rep(0, 55), 
    for(i in rownames(model_p)) {
      model_p[i,"Group_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] == x),i], na.rm = TRUE)
      model_p[i,"NotGroup_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] != x),i], na.rm = TRUE)
      
      # model_construct <- as.formula(paste0(i," ~ ImmuneCluster"))
      if(group_var == "Smoking_Status") {
        model_construct <- as.formula(paste0(i," ~ ImmuneCluster"))
      } else {
        model_construct<- as.formula(paste0(i," ~ ImmuneCluster + Smoking_Status + Site"))
      }
      model <- lme(model_construct,
                   random = ~ 1|Sample,
                   data = expr_subset,
                   na.action = na.omit)
      tTab <- summary(model)$tTable
      rownames(tTab)[2] <- "ImmuneCluster"
      model_p[i,"Beta"] <- tTab["ImmuneCluster","Value"]
      model_p[i,"P"] <- tTab["ImmuneCluster","p-value"]
    }
    model_p$Q <- p.adjust(model_p$P, n = nrow(model_p), method = "fdr")
    model_p$LogQ <- -log10(model_p$Q)
    model_p$Significant <- ((model_p$Q < 0.05) & (abs(model_p$Beta) > beta))
    model_p$SignificantPositive <- ((model_p$Q < 0.05) & (model_p$Beta > beta))
    model_p <- model_p %>% arrange(desc(SignificantPositive), Q)
  })
  names(each_model_res) <- clusters
  each_model_res
}

################################################################################
########################## Function to compare VAM scores of modules between clusters
vam_model_logistic <- function(scores, Cell_Type, k) {
  scores_immune <- scores[colnames(sce),]
  scores_immune$CellType <- sce$immune_CellType
  scores_immune$Smoking_Status <- sce$Smoking_Status
  scores_immune$Sample <- sce$Sample
  scores_immune$KCluster <- sce$immune_celda_cluster
  scores_immune$PCGA02_site <- sce$PCGA02_site
  scores_immune$BasalGroupMixed <- sce$BasalGroupMixed
  
  # scores_basal <- scores_immune[scores_immune$BasalGroupMixed %in% c("High Grade Basal Sample", "Low Grade Basal Sample"),]
  scores_basal <- scores_immune
  
  scores_basal <- scores_basal[scores_basal$CellType == Cell_Type,]
  scores_basal$CellType_vs_Not <- scores_basal$KCluster
  scores_basal$CellType_vs_Not <- unfactor(scores_basal$CellType_vs_Not)
  scores_basal$CellType_vs_Not[!scores_basal$CellType_vs_Not == k] <- 0
  scores_basal$CellType_vs_Not <- factor(scores_basal$CellType_vs_Not, levels = c(0, k), ordered = FALSE)
  # scores_basal <- within(scores_basal, CellType_vs_Not <- relevel(CellType_vs_Not, ref = 0))
  
  # scores_basal$BasalGroupMixed <- factor(scores_basal$BasalGroupMixed, levels = c("Low Grade Basal Sample", "High Grade Basal Sample"), ordered = FALSE)
  scores_basal$Smoking_Status <- factor(scores_basal$Smoking_Status, levels = levels(sce$Smoking_Status), ordered = FALSE)
  scores_basal$PCGA02_site <- factor(scores_basal$PCGA02_site, levels = unique(sce$PCGA02_site), ordered = FALSE)
  
  # scores_basal <- within(scores_basal, BasalGroupMixed <- relevel(BasalGroupMixed, ref = "Low Grade Basal Sample"))
  scores_basal <- within(scores_basal, Smoking_Status <- relevel(Smoking_Status, ref = "Former"))
  scores_basal <- within(scores_basal, PCGA02_site <- relevel(PCGA02_site, ref = "Roswell"))
  
  model_p <- data.frame(Module = colnames(scores),
                        KNot_Mean = rep(0, length(colnames(scores))),
                        K_Mean = rep(0, length(colnames(scores))),
                        # LGB_Mean = rep(0, 28),
                        # HGB_Mean = rep(0, 28),
                        Score_Logit = rep(0, length(colnames(scores))),
                        Score_P = rep(0, length(colnames(scores))),
                        Smoke_Logit = rep(length(colnames(scores))),
                        Smoke_P = rep(0, length(colnames(scores))),
                        row.names = colnames(scores))
  for(i in model_p$Module) {
    model_p[i, "KNot_Mean"] <- mean(scores_basal[which(scores_basal$CellType_vs_Not == 0),i])
    model_p[i, "K_Mean"] <- mean(scores_basal[which(scores_basal$CellType_vs_Not == k),i])
    
    # model_p[i, "LGB_Mean"] <- mean(scores_basal[which(scores_basal$BasalGroupMixed == "Low Grade Basal Sample"),i])
    # model_p[i, "HGB_Mean"] <- mean(scores_basal[which(scores_basal$BasalGroupMixed == "High Grade Basal Sample"),i])
    model_terms <- paste0("CellType_vs_Not ~ ", i, " + Smoking_Status + PCGA02_site + (1|Sample)")
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

################################################################################
########################### 4D: Macrophages

################### Look at DE external signatures
mac_modules_scores <- scores[,c(grep("Mac", colnames(scores), ignore.case = TRUE), grep("Myeloid", colnames(scores), ignore.case = TRUE))]
glmm_module_res_mac<- vam_model_logistic(scores = mac_modules_scores, Cell_Type = "Macrophages", k = 26)
mac_model_results <- glmm_module_res_mac$model_results %>% arrange(Score_Q)

saveRDS(mac_model_results, "./Objects_Results/Fig4D_DE_ExternalSignatures.rds")

################### Find modules
mac_cellprobs <- cellprobs[which(cellprobs$Cluster %in% c(26, 27)),]
diff_modules_mac <- immune_clust_de(expr = mac_cellprobs, clusters = unique(mac_cellprobs$Cluster), beta = 0.5)
diff_modules_res_mac <- Reduce(diff_modules_mac, f = rbind)
diff_modules_res_significant_mac <- diff_modules_res_mac[which(diff_modules_res_mac$SignificantPositive),]

saveRDS(diff_modules_res_significant_mac, "./Objects_Results/Fig4D_DE_Modules.rds")

################### Make heatmap of DE modules
sce_modular_mac <- sce_modular[,which(sce_modular$immune_CellType == "Macrophages")]

sce_modular_mac$immune_celda_cluster <- factor(sce_modular_mac$immune_celda_cluster, levels = c(26, 27))
cd <- as.data.frame(colData(sce_modular_mac))
# sce_modular_mac$Zilionis_hMac6 <- scale(sce_modular_mac$Zilionis_hMac6)

########## Subset to appropriate modules & reorder into desired order
lgcts <- assay(sce_modular_mac, "LogNormalize")
lgcts <- lgcts[unique(diff_modules_res_significant_mac$Module),]
lgcts <- lgcts[c( "L19", "L23", "L10", "L25", "L26"),]
lgcts <- t(scale(t(lgcts)))

########## Build heatmap annotations
column_ha_all <- HeatmapAnnotation(Celda_Cluster = sce_modular_mac$immune_celda_cluster,
                                   Basal_Group = sce_modular_mac$BasalGroupMixed,
                                   HMac6_Macrophages = sce_modular_mac$Zilionis_hMac6,
                                   col = list(Basal_Group = c("Low Grade Basal Sample" = "blue", "Mixed Grade Basal Sample" = "#E0B0FF", "High Grade Basal Sample" = "red",
                                                              "Too Few Basal Cells" = "#F0FFFF", "No Epithelial Cells" = "#F0FFFF"),
                                              Celda_Cluster = c("26" = "#bae8ca", "27" = "#17853d"),
                                              HMac6_Macrophages = circlize::colorRamp2(c(min(sce_modular_mac$Zilionis_hMac6), median(sce_modular_mac$Zilionis_hMac6), max(sce_modular_mac$Zilionis_hMac6)), c("#ffffff", "#f7cdd0", "#ce1256"))),
                                   annotation_legend_param = list(
                                     HMac6_Macrophages = list(
                                       title = "HMac6 Macrophages\nSignature VAM Score",
                                       at = c(min(sce_modular_mac$Zilionis_hMac6), median(sce_modular_mac$Zilionis_hMac6), max(sce_modular_mac$Zilionis_hMac6)),  # match breakpoints in the color function
                                       labels = c("Min", "Median", "Max")
                                     )),
                                   show_annotation_name = TRUE,
                                   show_legend = TRUE,
                                   annotation_name_rot = 0)


row_ha <- rowAnnotation(Cluster = diff_modules_res_significant_mac$Cluster,
                        col = list(Cluster = c("26" = "#bae8ca", "27" = "#17853d")),
                        show_legend = FALSE)

########## Adjust min/max for scaling of colors
lgcts[which(lgcts > 2)] <- 2
lgcts[which(lgcts < -2)] <- -2

################### Create heatmap
pdf(paste0("./Figures/Fig4D_Macrophage_Heatmap.pdf"), width = 9)
Heatmap(lgcts, 
        top_annotation = column_ha_all,
        right_annotation = row_ha,
        show_column_names = FALSE, 
        show_row_names = TRUE,
        column_split = sce_modular_mac$immune_celda_cluster,
        show_heatmap_legend = TRUE,
        # column_title = NULL,
        column_title_rot = 45,
        column_title_gp = gpar(fontsize = 10),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE)
dev.off()


################################################################################
########################### 4E: Neutrophils

################### Look at DE external signatures
neut_modules_scores <- scores[,c(grep("Neut", colnames(scores), ignore.case = TRUE), grep("NAN", colnames(scores), ignore.case = FALSE), 
                                 grep("TAN", colnames(scores), ignore.case = FALSE), grep("Xie", colnames(scores), ignore.case = FALSE))]

glmm_module_res_neut_k3 <- vam_model_logistic(scores = neut_modules_scores, Cell_Type = "Neutrophils", k = 3)
neut_model_results_k3 <- glmm_module_res_neut_k3$model_results %>% arrange(Score_Q)
saveRDS(neut_model_results_k3, "./Objects_Results/Fig4E_DE_ExternalSignatures_K3_vs_K1K4.rds")

glmm_module_res_neut_k1 <- vam_model_logistic(scores = neut_modules_scores, Cell_Type = "Neutrophils", k = 1)
neut_model_results_k1 <- glmm_module_res_neut_k1$model_results %>% arrange(Score_Q)
saveRDS(neut_model_results_k1, "./Objects_Results/Fig4E_DE_ExternalSignatures_K1_vs_K3K4.rds")

glmm_module_res_neut_k4 <- vam_model_logistic(scores = neut_modules_scores, Cell_Type = "Neutrophils", k = 4)
neut_model_results_k4 <- glmm_module_res_neut_k4$model_results %>% arrange(Score_Q)
saveRDS(neut_model_results_k4, "./Objects_Results/Fig4E_DE_ExternalSignatures_K4_vs_K1K3.rds")

# saveRDS(neut_model_results_k3, "./Objects_Results/Fig4E_DE_ExternalSignatures.rds")

######## Comparing scores between individual clusters
sce$immune_CellType_orig <- sce$immune_CellType
sce$immune_CellType[which(sce$immune_celda_cluster == 1)] <- "Other"
glmm_module_res_neut_k3_vs_k4 <- vam_model_logistic(scores = neut_modules_scores, Cell_Type = "Neutrophils", k = 3)
neut_model_results_k3_vs_k4 <- glmm_module_res_neut_k3_vs_k4$model_results %>% arrange(Score_Q)
saveRDS(neut_model_results_k3_vs_k4, "./Objects_Results/Fig4E_DE_ExternalSignatures_K3_vs_K4.rds")

sce$immune_CellType <- sce$immune_CellType_orig
sce$immune_CellType[which(sce$immune_celda_cluster == 3)] <- "Other"
glmm_module_res_neut_k1_vs_k4 <- vam_model_logistic(scores = neut_modules_scores, Cell_Type = "Neutrophils", k = 1)
neut_model_results_k1_vs_k4 <- glmm_module_res_neut_k1_vs_k4$model_results %>% arrange(Score_Q)
saveRDS(neut_model_results_k1_vs_k4, "./Objects_Results/Fig4E_DE_ExternalSignatures_K1_vs_K4.rds")

sce$immune_CellType <- sce$immune_CellType_orig
sce$immune_CellType[which(sce$immune_celda_cluster == 4)] <- "Other"
glmm_module_res_neut_k1_vs_k3 <- vam_model_logistic(scores = neut_modules_scores, Cell_Type = "Neutrophils", k = 1)
neut_model_results_k1_vs_k3 <- glmm_module_res_neut_k1_vs_k3$model_results %>% arrange(Score_Q)
saveRDS(neut_model_results_k1_vs_k3, "./Objects_Results/Fig4E_DE_ExternalSignatures_K1_vs_K3.rds")

################### Find modules
neut_cellprobs <- cellprobs[which(cellprobs$Cluster %in% c(1,3,4)),]
diff_modules_neut <- immune_clust_de(expr = neut_cellprobs, clusters = unique(neut_cellprobs$Cluster), beta = 0.35)
diff_modules_res_neut <- Reduce(diff_modules_neut, f = rbind)
diff_modules_res_significant_neut <- diff_modules_res_neut[which(diff_modules_res_neut$SignificantPositive),]

saveRDS(diff_modules_res_significant_neut, "./Objects_Results/Fig4E_DE_Modules.rds")

######## Comparing individual clusters 
K1vsK4_cellprobs <- cellprobs[which(cellprobs$Cluster %in% c(1,4)),]
diff_modules_K1vsK4 <- immune_clust_de(expr = K1vsK4_cellprobs, clusters = unique(K1vsK4_cellprobs$Cluster), beta = 0.35)
diff_modules_res_K1vsK4 <- Reduce(diff_modules_K1vsK4, f = rbind)
diff_modules_res_significant_K1vsK4 <- diff_modules_res_K1vsK4[which(diff_modules_res_K1vsK4$SignificantPositive),]
saveRDS(diff_modules_res_significant_K1vsK4, "./Objects_Results/Fig4E_DE_K1vsK4_Modules.rds")

K1vsK3_cellprobs <- cellprobs[which(cellprobs$Cluster %in% c(1,3)),]
diff_modules_K1vsK3 <- immune_clust_de(expr = K1vsK3_cellprobs, clusters = unique(K1vsK3_cellprobs$Cluster), beta = 0.35)
diff_modules_res_K1vsK3 <- Reduce(diff_modules_K1vsK3, f = rbind)
diff_modules_res_significant_K1vsK3 <- diff_modules_res_K1vsK3[which(diff_modules_res_K1vsK3$SignificantPositive),]
saveRDS(diff_modules_res_significant_K1vsK3, "./Objects_Results/Fig4E_DE_K1vsK3_Modules.rds")

K3vsK4_cellprobs <- cellprobs[which(cellprobs$Cluster %in% c(3,4)),]
diff_modules_K3vsK4 <- immune_clust_de(expr = K3vsK4_cellprobs, clusters = unique(K3vsK4_cellprobs$Cluster), beta = 0.35)
diff_modules_res_K3vsK4 <- Reduce(diff_modules_K3vsK4, f = rbind)
diff_modules_res_significant_K3vsK4 <- diff_modules_res_K3vsK4[which(diff_modules_res_K3vsK4$SignificantPositive),]
saveRDS(diff_modules_res_significant_K3vsK4, "./Objects_Results/Fig4E_DE_K3vsK4_Modules.rds")

################### Make heatmap of DE modules
sce_modular_neut <- sce_modular[,which(sce_modular$immune_CellType == "Neutrophils")]

sce_modular_neut$immune_celda_cluster <- factor(sce_modular_neut$immune_celda_cluster, levels = c(1, 3, 4))
cd <- as.data.frame(colData(sce_modular_neut))

########## Subset to appropriate modules & reorder into desired order
lgcts <- assay(sce_modular_neut, "LogNormalize")
# lgcts <- lgcts[unique(diff_modules_res_significant_neut$Module),]
# lgcts <- lgcts[c("L1", "L10", "L25", "L26", "L28",
#                  "L7", "L15", "L13",
#                  "L34", "L36"),]
lgcts <- lgcts[unique(c(Fig4E_DE_K1vsK4_Modules$Module, Fig4E_DE_K1vsK3_Modules$Module, Fig4E_DE_K3vsK4_Modules$Module)),]
lgcts <- lgcts[c("L7", "L37", "L10", "L25", "L26", 
                 "L1", "L28", 
                 "L8", "L43", 
                 "L11", "L12", "L13", "L15", 
                 "L14", "L34", "L36"
                 ),]
lgcts <- t(scale(t(lgcts)))

########## Build heatmap annotations
column_ha_all <- HeatmapAnnotation(Celda_Cluster = sce_modular_neut$immune_celda_cluster,
                                   Basal_Group = sce_modular_neut$BasalGroupMixed,
                                   Salcher_NAN = sce_modular_neut$NAN_5B,
                                   Salcher_TAN = sce_modular_neut$TAN_5B,
                                   hNeutro5 = sce_modular_neut$Zilionis_hNeutro5,
                                   col = list(Basal_Group = c("Low Grade Basal Sample" = "blue", "Mixed Grade Basal Sample" = "#E0B0FF", "High Grade Basal Sample" = "red",
                                                              "Too Few Basal Cells" = "#F0FFFF", "No Epithelial Cells" = "#F0FFFF"),
                                              Celda_Cluster = c("1" = "lightpink", "3" = "#DE3163", "4" = "#9F2B68"),
                                              Salcher_NAN = circlize::colorRamp2(c(min(sce_modular_neut$NAN_6B), median(sce_modular_neut$NAN_6B), max(sce_modular_neut$NAN_6B)), c("#ffffff", "#fee8c8", "#e34a33")),
                                              Salcher_TAN = circlize::colorRamp2(c(min(sce_modular_neut$TAN_6B), median(sce_modular_neut$TAN_6B), max(sce_modular_neut$TAN_6B)), c("#ffffff", "#a6dba0", "#1b7837")),
                                              hNeutro5 = circlize::colorRamp2(c(min(sce_modular_neut$Zilionis_hNeutro5), median(sce_modular_neut$Zilionis_hNeutro5), max(sce_modular_neut$Zilionis_hNeutro5)), c("#ffffff", "#bdd7e7", "#08519c"))), # "2" = "#FF69B4", 
                                   annotation_legend_param = list(
                                     Salcher_NAN = list(
                                       title = "Salcher_NAN VAM Score",
                                       at = c(min(sce_modular_neut$NAN_5B), median(sce_modular_neut$NAN_5B), max(sce_modular_neut$NAN_5B)),  # match breakpoints in the color function
                                       labels = c("Min", "Median", "Max")
                                     ),
                                     Salcher_TAN = list(
                                       title = "Salcher_TAN VAM Score",
                                       at = c(min(sce_modular_neut$TAN_5B), median(sce_modular_neut$TAN_5B), max(sce_modular_neut$TAN_5B)),  # match breakpoints in the color function
                                       labels = c("Min", "Median", "Max")
                                     ),
                                     hNeutro5 = list(
                                       title = "Zilionis hNeutro5 VAM Score",
                                       at = c(min(sce_modular_neut$Zilionis_hNeutro5), median(sce_modular_neut$Zilionis_hNeutro5), max(sce_modular_neut$Zilionis_hNeutro5)),  # match breakpoints in the color function
                                       labels = c("Min", "Median", "Max")
                                     )
                                   ),
                                   show_annotation_name = TRUE,
                                   show_legend = TRUE,
                                   annotation_name_rot = 0)

########## Create a binary matrix for row annotations
num_modules <- length(unique(rownames(lgcts)))
num_clusters <- 3
annotation_matrix <- matrix(0, nrow = num_modules, ncol = num_clusters)

rownames(annotation_matrix) <- unique(rownames(lgcts))
# colnames(annotation_matrix) <- c(1,3,4)
# 
# # Fill in the matrix based on enrichment
# for (i in 1:num_modules) {
#   # Get the current module name
#   module_num <- rownames(annotation_matrix)[i]
#   
#   # Find rows corresponding to the current module
#   module_rows <- diff_modules_res_significant_neut[diff_modules_res_significant_neut$Module == module_num, ]
#   
#   # Extract the clusters as numeric values
#   enriched_clusters <- as.numeric(as.character(module_rows$Cluster))
#   
#   # Filter out NA values
#   enriched_clusters <- enriched_clusters[!is.na(enriched_clusters)]
#   
#   # Update the annotation matrix if enriched clusters exist
#   if (length(enriched_clusters) > 0) {
#     annotation_matrix[i, which(colnames(annotation_matrix) %in% enriched_clusters)] <- 1
#   }
# }

########## Convert to a data frame for better handling
annotation_df <- as.data.frame(annotation_matrix)
# colnames(annotation_df) <- paste0("Cluster_", colnames(annotation_df))

# annotation_df[,4] <- 0
# annotation_df[,5] <- 0
# annotation_df[,6] <- 0

colnames(annotation_df)[1] <- "K1_vs_K3"
colnames(annotation_df)[2] <- "K1_vs_K4"
colnames(annotation_df)[3] <- "K3_vs_K4"

annotation_df[which(rownames(annotation_df) %in% Fig4E_DE_K1vsK3_Modules$Module[which(Fig4E_DE_K1vsK3_Modules$Cluster == 1)]),1] <- 1
annotation_df[which(rownames(annotation_df) %in% Fig4E_DE_K1vsK3_Modules$Module[which(Fig4E_DE_K1vsK3_Modules$Cluster == 3)]),1] <- 2

annotation_df[which(rownames(annotation_df) %in% Fig4E_DE_K1vsK4_Modules$Module[which(Fig4E_DE_K1vsK4_Modules$Cluster == 1)]),2] <- 1
annotation_df[which(rownames(annotation_df) %in% Fig4E_DE_K1vsK4_Modules$Module[which(Fig4E_DE_K1vsK4_Modules$Cluster == 4)]),2] <- 2

annotation_df[which(rownames(annotation_df) %in% Fig4E_DE_K3vsK4_Modules$Module[which(Fig4E_DE_K3vsK4_Modules$Cluster == 3)]),3] <- 1
annotation_df[which(rownames(annotation_df) %in% Fig4E_DE_K3vsK4_Modules$Module[which(Fig4E_DE_K3vsK4_Modules$Cluster == 4)]),3] <- 2

row_ha <- rowAnnotation(
  # K1 = annotation_df$Cluster_1,
  # K2 = annotation_df$Cluster_2,
  # K3 = annotation_df$Cluster_3,
  # K4 = annotation_df$Cluster_4,
  K1_vs_K3 = annotation_df$K1_vs_K3,
  K1_vs_K4 = annotation_df$K1_vs_K4,
  K3_vs_K4 = annotation_df$K3_vs_K4,
  col = list(
    # K1 = c("1" = "lightpink", "0" = "white"),
    # K2 = c("1" = "#FF69B4", "0" = "white"),
    # K3 = c("1" = "#DE3163", "0" = "white"),
    # K4 = c("1" = "#9F2B68", "0" = "white"),
    K1_vs_K3 = c("1" = "#B39DDB", "2" = "#673AB7", "0" = "white"),
    K1_vs_K4 = c("1" = "#81D4FA", "2" = "#0288D1", "0" = "white"),
    K3_vs_K4 = c("1" = "#80CBC4", "2" = "#00695C", "0" = "white")
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "bottom",
  show_legend = TRUE # Suppress legend
)

########## Adjust min/max for scaling of colors
lgcts[which(lgcts > 2)] <- 2
lgcts[which(lgcts < -2)] <- -2

################### Create heatmap
pdf(paste0("./Figures/Fig4E_Neutrophil_Heatmap_temp.pdf"), width = 9)
Heatmap(lgcts, 
        top_annotation = column_ha_all,
        right_annotation = row_ha,
        show_column_names = FALSE, 
        show_row_names = TRUE,
        column_split = sce_modular_neut$immune_celda_cluster,
        show_heatmap_legend = TRUE,
        # column_title = NULL,
        column_title_rot = 45,
        column_title_gp = gpar(fontsize = 10),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE)
dev.off()


################################################################################
###########################  4F: B Cells

################### Look at DE external signatures
bcell_modules_scores <- scores[,c(grep("Hao", colnames(scores), ignore.case = TRUE), grep("Madissoon_B", colnames(scores), ignore.case = TRUE),
                                  grep("BCell", colnames(scores), ignore.case = TRUE), grep("B_Cell", colnames(scores), ignore.case = TRUE))]
glmm_module_res_bcell <- vam_model_logistic(scores = bcell_modules_scores, Cell_Type = "B Cells", k = 13)
bcell_model_results <- glmm_module_res_bcell$model_results %>% arrange(Score_Q)

saveRDS(bcell_model_results, "./Objects_Results/Fig4F_DE_ExternalSignatures.rds")

################### Find modules
b_cellprobs <- cellprobs[which(cellprobs$Cluster %in% c(7, 13)),]
diff_modules_bcell <- immune_clust_de(expr = b_cellprobs, clusters = unique(b_cellprobs$Cluster), beta = 0.5)
diff_modules_res_bcell <- Reduce(diff_modules_bcell, f = rbind)
diff_modules_res_significant_bcell <- diff_modules_res_bcell[which(diff_modules_res_bcell$SignificantPositive),]

saveRDS(diff_modules_res_significant_bcell, "./Objects_Results/Fig4F_DE_Modules.rds")

################### Find DE IG genes
sce_ig <- sce@int_colData@listData[["altExps"]]@listData[["IG"]]@se
sce_ig$immune_celda_cluster <- sce$immune_celda_cluster

sce_ig <- runDEAnalysis(inSCE = sce_ig, method = "wilcox", useAssay = "logcounts", 
                        class = "immune_celda_cluster", classGroup1 = c(7), classGroup2 = c(13), 
                        groupName1 = "c7", groupName2 = "c13", analysisName = "c7_VS_c13", overwrite = TRUE)
DEG <- getDEGTopTable(sce_ig, "c7_VS_c13")
head(DEG)

################### Make heatmap of DE modules
sce_modular_bcell <- sce_modular[,which(sce_modular$immune_CellType == "B Cells")]

sce_modular_bcell$immune_celda_cluster <- factor(sce_modular_bcell$immune_celda_cluster, levels = c(7, 13))
cd <- as.data.frame(colData(sce_modular_bcell))

########## Subset to appropriate modules & reorder into desired order
lgcts <- assay(sce_modular_bcell, "LogNormalize")
lgcts <- lgcts[unique(diff_modules_res_significant_bcell$Module),]
lgcts <- lgcts[c("L34", "L2", "L3", "L7", "L9", "L16", "L18", "L20", "L21", "L22", "L24", "L25", "L53", "L56", "L57"),]
lgcts <- t(scale(t(lgcts)))

########## Build heatmap annotations
column_ha_all <- HeatmapAnnotation(Celda_Cluster = sce_modular_bcell$immune_celda_cluster,
                                   Basal_Group = sce_modular_bcell$BasalGroupMixed,
                                   B_Cell_Activation = sce_modular_bcell$Cell_B_cell_activation,
                                   col = list(Basal_Group = c("Low Grade Basal Sample" = "blue", "Mixed Grade Basal Sample" = "#E0B0FF", "High Grade Basal Sample" = "red",
                                                              "Too Few Basal Cells" = "#F0FFFF", "No Epithelial Cells" = "#F0FFFF"),
                                              Celda_Cluster = c("7" = "#f0c189", "13" = "#c26d08"),
                                              B_Cell_Activation = circlize::colorRamp2(c(min(sce_modular_bcell$Cell_B_cell_activation), median(sce_modular_bcell$Cell_B_cell_activation), max(sce_modular_bcell$Cell_B_cell_activation)), c("#ffffff", "#a6bddb", "#2b8cbe"))),
                                   annotation_legend_param = list(
                                     B_Cell_Activation = list(
                                       title = "B Cell Ativation\nSignature VAM Score",
                                       at = c(min(sce_modular_bcell$Cell_B_cell_activation), median(sce_modular_bcell$Cell_B_cell_activation), max(sce_modular_bcell$Cell_B_cell_activation)),  # match breakpoints in the color function
                                       labels = c("Min", "Median", "Max")
                                     )),
                                   show_annotation_name = TRUE,
                                   show_legend = TRUE,
                                   annotation_name_rot = 0)


row_ha <- rowAnnotation(Cluster = c(7, rep(13, 14)), #diff_modules_res_significant_bcell$Cluster,
                        col = list(Cluster = c("7" = "#f0c189", "13" = "#c26d08")),
                        show_legend = FALSE)

##########  Adjust min/max for scaling of colors
lgcts[which(lgcts > 2)] <- 2
lgcts[which(lgcts < -2)] <- -2

################### Create heatmap
pdf(paste0("./Figures/Fig4F_BCell_Heatmap.pdf"), width = 9)
Heatmap(lgcts, 
        top_annotation = column_ha_all,
        right_annotation = row_ha,
        show_column_names = FALSE, 
        show_row_names = TRUE,
        column_split = sce_modular_bcell$immune_celda_cluster,
        show_heatmap_legend = TRUE,
        # column_title = NULL,
        column_title_rot = 45,
        column_title_gp = gpar(fontsize = 10),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE)
dev.off()

