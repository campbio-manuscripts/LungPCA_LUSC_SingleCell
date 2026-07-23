library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(singleCellTK)
library(celda)
library(ggplot2)
library(nlme)
library(emmeans)
library(pheatmap)
library(ComplexHeatmap)
library(ggpubr)
library(scater)
library(stringr)
library(GSVA)
library(biomaRt)
library(limma)
library(hypeR)
library(reshape2)
library(edgeR)
library(sccomp)
library(scales)
library(gridExtra)
library(xlsx)

################################################################################
######### Read in SCE
sce <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/Fig4-6_Paper_Scripts_Results/Misc_Scripts/EpithelialBrushBiopsyCells_SCE_Final.rds")
sce$epithelial_CellType <- gsub("^Basal Cells", "LG Basal Cells", sce$epithelial_CellType)
sce_orig <- sce

celltype_colors <- readRDS("../celltype_colors.rds")
smoking_colors <- readRDS("../smoking_colors.rds")
histology_colors <- readRDS("../histology_colors.rds")
sampletype_cols <- readRDS("../sampletype_cols.rds")
names(sampletype_cols)[1] <- "Biopsy"

##############################################################################################################
##################################### Figure 5A
pdf("./Figures/Fig5A_CellType_UMAP.pdf", width = 6, height = 6)
plotSCEDimReduceColData(sce, 
                        colorBy = "epithelial_CellType", 
                        reducedDimName = "celda_UMAP", 
                        colorScale = celltype_colors,
                        dotSize = .7,
                        labelClusters = FALSE) + 
  theme_classic() + 
  guides(color="none") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) +
  xlab("UMAP1") + 
  ylab("UMAP2")
dev.off()

pdf("./Figures/Fig5A_Smoking_UMAP_Final.pdf", width = 6, height = 6)
plotSCEDimReduceColData(sce, 
                        colorBy = "Smoking_Status", 
                        reducedDimName = "celda_UMAP", 
                        colorScale = smoking_colors,
                        dotSize = .7,
                        labelClusters = FALSE) + 
  theme_classic() + 
  guides(color="none") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) +
  xlab("UMAP1") + 
  ylab("UMAP2")
dev.off()

pdf("./Figures/Fig5A_SampleType_UMAP.pdf", width = 6, height = 6)
plotSCEDimReduceColData(sce, 
                        colorBy = "SampleType", 
                        reducedDimName = "celda_UMAP", 
                        colorScale = sampletype_cols,
                        dotSize = .7,
                        labelClusters = FALSE) + 
  theme_classic() + 
  guides(color="none") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) +
  xlab("UMAP1") + 
  ylab("UMAP2")
dev.off()


##############################################################################################################
##################################### Figure 5B

sce <- sce_orig

################### Organize object
f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")), 
                                    colData = colData(sce))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")

cellprobs <- data.frame(t(assay(sce_modular,"LogNormalize")))
cellprobs$Cluster <- sce_modular$epithelial_CellType 
mean_cellprobs <- cellprobs %>% group_by(Cluster) %>% summarize_all(mean)
mean_cellprobs <- data.frame(t(mean_cellprobs))
colnames(mean_cellprobs) <- mean_cellprobs["Cluster",]
mean_cellprobs <- mean_cellprobs[2:71,] ## 
mean_cellprobs <- sapply(mean_cellprobs, as.numeric)
mean_cellprobs <- data.frame(mean_cellprobs, row.names = colnames(cellprobs)[1:70])
cellprobs$Sample <- sce$Sample
cellprobs$Smoking_Status <- sce$Smoking_Status

################### Find modules
epithelial_clust_de <- function(expr = cellprobs, group_var = "Cluster", clusters) {
  expr_subset <- expr
  each_model_res <- lapply(clusters, function(x) {
    expr_subset[,"epithelialCluster"] <- ifelse(expr_subset[,group_var] == x, "B", "A")
    model_p <- data.frame(Cluster = x, Module = names(expr_subset)[1:70], P = rep(0, 70), Beta = rep(0, 70), Group_Mean = rep(0, 70), NotGroup_Mean = rep(0, 70), row.names = names(expr_subset)[1:70]) # Smoking_P = rep(0, 55), Smoking_Beta = rep(0, 55), 
    for(i in rownames(model_p)) {
      model_p[i,"Group_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] == x),i], na.rm = TRUE)
      model_p[i,"NotGroup_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] != x),i], na.rm = TRUE)
      
      # model_construct <- as.formula(paste0(i," ~ epithelialCluster"))
      if(group_var == "Smoking_Status") {
        model_construct <- as.formula(paste0(i," ~ epithelialCluster"))
      } else {
        model_construct<- as.formula(paste0(i," ~ epithelialCluster + Smoking_Status"))
      }
      model <- lme(model_construct,
                   random = ~ 1|Sample,
                   data = expr_subset,
                   control = lmeControl(opt = "optim"),
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

diff_modules <- epithelial_clust_de(clusters = unique(cellprobs$Cluster))
diff_modules_res <- Reduce(diff_modules, f = rbind)
diff_modules_res_significant <- diff_modules_res[which(diff_modules_res$SignificantPositive),]

saveRDS(diff_modules_res_significant, "./Objects_Results/Fig5B_diffExpModules_byCellType_Significant_Final.rds")

################### Specify order of cell types
celltype_order <- c("Ciliated Cells", "LG Basal Cells", "HG Basal Cells", "Goblet Cells", "Perigoblet Cells",
                    "Club Cells", "Nasal Secretory Cells", "Squamous Hillock Cells", "Mucous SMG Cells", "Serous SMG Cells", "Keratinizing Epithelial Cells")
celltype_modules <- as.data.frame(matrix(NA, ncol = 2, nrow = length(unique(sce$epithelial_CellType))))
celltype_modules$V1 <- celltype_order

for (i in 1:nrow(celltype_modules)) {
  celltype <- diff_modules_res_significant[which(diff_modules_res_significant$Cluster == celltype_modules$V1[i]),]
  celltype_modules$V2[i] <- celltype$Module[which(celltype$Beta == max(celltype$Beta))]
}

################### Prep heatmap data
m <- celltype_modules$V2
m <- gsub("L66", "L62", m) ## choose a different HGB module
m <- gsub("L31", "L29", m) ## choose a different PG module

scaled_expr <- t(assay(sce_modular,"LogNormalize")[m,])
scaled_expr <- scale(scaled_expr)

colnames(scaled_expr)[which(colnames(scaled_expr) == "L5")] <- "L5 (TPPP3)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L70")] <- "L70 (KRT5)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L62")] <- "L62 (KRT6A)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L15")] <- "L15 (MUC5AC)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L29")] <- "L29 (TXN, SERPINB1)"
# colnames(scaled_expr)[which(colnames(scaled_expr) == "L22")] <- "L22 (CEACAM5)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L10")] <- "L10 (SCGB3A1)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L25")] <- "L25 (LYPD2)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L60")] <- "L60 (KRT13)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L54")] <- "L54 (CST1)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L14")] <- "L14 (MSMB)"
colnames(scaled_expr)[which(colnames(scaled_expr) == "L59")] <- "L59 (S100A9)"

df <- cbind(data.frame(CellType = sce_modular$epithelial_celda_cluster), scaled_expr)
df_markersummary <- df %>% group_by(CellType) %>% summarize_all(mean)
df_markersummary_pivot <- df_markersummary %>% pivot_longer(!CellType, names_to = "Module", values_to = "Expression")
df_markersummary_pivot$CellType <- paste0("K", df_markersummary_pivot$CellType)
df_markersummary_pivot$CellType <- factor(df_markersummary_pivot$CellType,
                                          levels = c("K6", "K20", "K31", "K28", "K35",
                                                     "K37", "K30", "K27", "K34", "K36",
                                                     "K21", "K38", "K23", "K40", "K22", 
                                                     "K39", "K29", "K32", "K33", "K24",
                                                     "K26", "K25", "K16", "K7", "K18",
                                                     "K15", "K17", "K10", "K13", "K19",
                                                     "K14", "K9", "K12", "K8", "K11", "K4",
                                                     "K2", "K3", "K1", "K5"),
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

pdf("./Figures/Fig5B_CellType_Modules_Heatmap_Final.pdf", width = 4, height = 7)
marker_scalehm
dev.off()


##############################################################################################################
##################################### Figure 5C
sce <- sce_orig

sce$SampleType <- factor(sce$SampleType, levels = c("Biopsy", "BronchialBrush", "NasalBrush"))
sce$epithelial_celda_cluster <- paste0("K", sce$epithelial_celda_cluster, sep = "")
sce$Smoking_Status <- factor(sce$Smoking_Status, levels = c("Former", "Current"), ordered = FALSE)
sce <- sce[,-which(is.na(sce$Smoking_Status))]
sce$PCGA02_site <- factor(sce$PCGA02_site, levels = c("Roswell", "UCL"))
sce$Sample <- paste0("S", sce$Sample, sep = "")

################### Run SCComp on celda cluster
# model_with_factor_association = 
#   sce |>
#   sccomp_estimate(  #sccomp_glm
#     formula_composition = ~ 1 + SampleType + Smoking_Status + PCGA02_site,
#     # formula_variability = ~ 1,
#     .sample = Sample, 
#     .cell_group = epithelial_celda_cluster, 
#     # bimodal_mean_variability_association = TRUE,
#     # cores = 1, 
#     enable_loo = FALSE,
#     inference_method = "hmc"
#     # variational_inference = FALSE
#   ) |> 
#   sccomp_remove_outliers(cores = 4) |> 
#   sccomp_test()

estimate_raw <- sce |>
  sccomp_estimate(
    formula_composition = ~ 1 + SampleType + Smoking_Status + PCGA02_site,
    .sample = Sample,
    .cell_group = epithelial_celda_cluster,
    enable_loo = FALSE,
    inference_method = "hmc",
    cores = 4,
    verbose = TRUE
  )

model_with_factor_association <- estimate_raw |>
  sccomp_test()

saveRDS(model_with_factor_association, "./Objects_Results/sccomp_results_SampleType_covSmokingSite_Final.rds")

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

stats_values$epithelial_celda_cluster <- factor(stats_values$epithelial_celda_cluster, levels = c("K6", "K20", "K31", "K28", "K35",
                                                                                                  "K37", "K30", "K27", "K34", "K36",
                                                                                                  "K21", "K38", "K23", "K40", "K22", 
                                                                                                  "K39", "K29", "K32", "K33", "K24",
                                                                                                  "K26", "K25", "K16", "K7", "K18",
                                                                                                  "K15", "K17", "K10", "K13", "K19",
                                                                                                  "K14", "K9", "K12", "K8", "K11", "K4",
                                                                                                  "K2", "K3", "K1", "K5"))

stats_values <- stats_values[-which(stats_values$parameter == "(Intercept)"),]
stats_values <- stats_values[-which(stats_values$parameter == "PCGA02_siteUCL"),]
stats_values$parameter <- gsub("SampleTypeBronchialBrush", "BBr", stats_values$parameter)
stats_values$parameter <- gsub("SampleTypeNasalBrush", "NBr", stats_values$parameter)
stats_values$parameter <- gsub("Smoking_StatusCurrent", "SmCurr", stats_values$parameter)

################### Plot
########## Bubbleplot
colors <- colorRampPalette(c("#053061", "#053061", "#053061", "#2166AC", "#4393C3", 
                             "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#67001F", 
                             "#67001F", "#67001F"))(200)

rescales <- rescale(c(min(na.omit(stats_values$c_effect)), 0, max(na.omit(stats_values$c_effect))), 
                    limits = c(min(na.omit(stats_values$c_effect)), max(na.omit(stats_values$c_effect))))

p <- ggplot(stats_values, aes(x = parameter, y = epithelial_celda_cluster, label = sig, color = c_effect, 
                              size = ifelse(is.na(logp), 2, logp), shape = ifelse(is.na(c_effect), "Missing", "Present"))) + 
  geom_point(pch=1, colour = "black", stroke = 2) +
  geom_point(stat = "identity", stroke = 2, na.rm = T, shape = 19) + # colour="black",  stroke = 2, 
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5)) +
  labs(x = "", y = "", color = "Effect Size") + 
  guides(size = F, shape = F) +
  geom_text(color = "white", size = 3, vjust = 0.9) +
  scale_shape_manual(values = c(Missing = 4, Present = 19)) +
  scale_color_gradientn(colours = colors, values = rescales) +
  scale_x_discrete(labels = label_wrap(10)) +
  guides(color="none")

########## Boxplots
cd <- as.data.frame(colData(sce))
cd$epithelial_celda_cluster <- factor(cd$epithelial_celda_cluster, levels = c("K6", "K20", "K31", "K28", "K35",
                                                                              "K37", "K30", "K27", "K34", "K36",
                                                                              "K21", "K38", "K23", "K40", "K22", 
                                                                              "K39", "K29", "K32", "K33", "K24",
                                                                              "K26", "K25", "K16", "K7", "K18",
                                                                              "K15", "K17", "K10", "K13", "K19",
                                                                              "K14", "K9", "K12", "K8", "K11", "K4",
                                                                              "K2", "K3", "K1", "K5"))

cd_sample <- cd %>% dplyr::group_by(epithelial_celda_cluster) %>% dplyr::count(SampleType) %>% na.omit()

sample <- ggplot(cd_sample, aes(y = epithelial_celda_cluster, fill = SampleType, x = n)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  scale_fill_manual(values = c("Biopsy" = "#0C7BDC", "BronchialBrush" = "#FFC20A", "NasalBrush" = "#00AB41")) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(fill = "none") +
  labs(x = "") +
  scale_x_continuous(labels = scales::percent)

cd_smoke <- cd %>% dplyr::group_by(epithelial_celda_cluster) %>% dplyr::count(Smoking_Status) %>% na.omit()

smoke <- ggplot(cd_smoke, aes(y = epithelial_celda_cluster, fill = Smoking_Status, x = n)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  scale_fill_manual(values = c("Former" = "grey40", "Current" = "red")) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(fill = "none") +
  labs(x = "") +
  scale_x_continuous(labels = scales::percent)


################### Combine & save plots
pdf(paste0("./Figures/Fig5C_bubbleplot_SampleType_covSmokingSite_separateSmokingBoxplots_Final.pdf"), height = 8, width = 6)
grid.arrange(p, sample, smoke, ncol = 3, nrow=1, widths = c(1.5, 1, 1))
dev.off()


##############################################################################################################
##################################### Figure 5D & E

sce <- sce_orig

########################## Prep data
f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")), 
                                    colData = colData(sce))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")

cellprobs <- data.frame(t(assay(sce_modular,"LogNormalize")))
cellprobs$Cluster <- sce_modular$epithelial_celda_cluster ## 
cellprobs$CellType <- sce_modular$epithelial_CellType
cellprobs$Sample <- sce_modular$Sample
cellprobs$Subject <- sce_modular$Subject
cellprobs$Smoking_Status <- sce_modular$Smoking_Status
cellprobs$SampleType <- sce_modular$SampleType
cellprobs$Site <- sce_modular$PCGA02_site

########################## Identify DE modules
epi_clust_de <- function(expr = cellprobs, group_var = "Cluster", clusters) {
  expr_subset <- expr
  each_model_res <- lapply(clusters, function(x) {
    expr_subset[,"EpiCluster"] <- ifelse(expr_subset[,group_var] == x, "B", "A")
    model_p <- data.frame(Cluster = x, Module = names(expr_subset)[1:70], P = rep(0, 70), Beta = rep(0, 70), Group_Mean = rep(0, 70), NotGroup_Mean = rep(0, 70), row.names = names(expr_subset)[1:70]) # Smoking_P = rep(0, 55), Smoking_Beta = rep(0, 55), 
    for(i in rownames(model_p)) {
      model_p[i,"Group_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] == x),i], na.rm = TRUE)
      model_p[i,"NotGroup_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] != x),i], na.rm = TRUE)
      
      if(group_var == "SampleType") {
        model_construct <- as.formula(paste0(i," ~ EpiCluster"))
      } else {
        model_construct<- as.formula(paste0(i," ~ EpiCluster + Smoking_Status + Site"))
      }
      model <- lme(model_construct,
                   random = ~ 1|Sample,
                   data = expr_subset,
                   control = lmeControl(opt = "optim"),
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


##################################### Figure 5D: Nasal secretory cell type vs others
# sec_cellprobs <- cellprobs[which(cellprobs$CellType %in% c("Club Cells", "Perigoblet Cells", "Goblet Cells", "Nasal Secretory Cells")),]
# diff_modules_sec <- epi_clust_de(expr = sec_cellprobs, group_var = "CellType", clusters = unique(sec_cellprobs$CellType))
# diff_modules_res_sec <- Reduce(diff_modules_sec, f = rbind)
# diff_modules_res_significant_sec <- diff_modules_res_sec[which(diff_modules_res_sec$SignificantPositive),]
# 
# saveRDS(diff_modules_res_sec, "./Objects_Results/Fig5D_DEModules_SecretoryCellTypes_Final.rds")
# saveRDS(diff_modules_res_significant_sec, "./Objects_Results/Fig5D_DEModules_SecretoryCellTypes_Significant_Final.rds")

diff_modules_res_sec_celltypes <- readRDS("./Objects_Results/Fig5D_DEModules_SecretoryCellTypes_Final.rds")
diff_modules_res_significant_sec_celltypes <- readRDS("./Objects_Results/Fig5D_DEModules_SecretoryCellTypes_Significant_Final.rds")

################### Make heatmap

########## Modules upregulated in Nasal Secretory Cells as a whole
modules <- unique(diff_modules_res_significant_sec_celltypes$Module[which(diff_modules_res_significant_sec_celltypes$Beta > 0.7)])
modules <- unique(diff_modules_res_significant_sec_celltypes$Module)

########## The factorizeMatrix function generates a "factorized" matrix between features/feature modules and samples/cell clusters/cells.
factoMat <- factorizeMatrix(sce, useAssay = "decontXcounts")
cellModuleMatrix <- factoMat$counts$cell

########## Standard normalization of counts
pop <- SingleCellExperiment(assays = SimpleList("factorizeCounts" = cellModuleMatrix), colData = colData(sce))
pop <- logNormCounts(pop, assay.type = "factorizeCounts", name = "logcounts")

pop <- pop[,which(pop$epithelial_celda_cluster %in% c(21:27, 29, 30, 32:40))]
pop$epithelial_celda_cluster <- factor(pop$epithelial_celda_cluster, levels = c(21:23, 38:40, 24:26, 29, 32, 33, 27, 30, 34:37))
pop$epithelial_CellType <- factor(pop$epithelial_CellType, levels = c("Club Cells", "Goblet Cells", "Perigoblet Cells", "Nasal Secretory Cells"))

########## LogFC threshold
logFC <- as.data.frame(matrix(NA, nrow = length(modules)), ncol = 5)
logFC[,1] <- modules
for (i in 1:length(modules)) {
  mod_club <- mean(logcounts(pop)[modules[i],which(pop$epithelial_CellType == "Club Cells")])
  mod_other <- mean(logcounts(pop)[modules[i],which(!pop$epithelial_CellType == "Club Cells")])
  logFC[i,2] <- mod_club - mod_other
  
  mod_goblet <- mean(logcounts(pop)[modules[i],which(pop$epithelial_CellType == "Goblet Cells")])
  mod_other <- mean(logcounts(pop)[modules[i],which(!pop$epithelial_CellType == "Goblet Cells")])
  logFC[i,3] <- mod_goblet - mod_other
  
  mod_pg <- mean(logcounts(pop)[modules[i],which(pop$epithelial_CellType == "Perigoblet Cells")])
  mod_other <- mean(logcounts(pop)[modules[i],which(!pop$epithelial_CellType == "Perigoblet Cells")])
  logFC[i,4] <- mod_pg - mod_other
  
  mod_nas <- mean(logcounts(pop)[modules[i],which(pop$epithelial_CellType == "Nasal Secretory Cells")])
  mod_other <- mean(logcounts(pop)[modules[i],which(!pop$epithelial_CellType == "Nasal Secretory Cells")])
  logFC[i,5] <- mod_nas - mod_other
}
colnames(logFC) <- c("Module", "Club", "Goblet", "Perigoblet", "Nasal Secretory")

########## Subset to modules that have LogFC > 1.25 in their designated cell type
modules <- c("L10", "L9", "L13", "L56", "L2",
             "L12", "L15", "L14", "L55", "L17",
             "L29", "L22",
             "L25", "L23", "L24", "L32", "L59")
# modules <- logFC$Module[which(abs(logFC$LogFC) > 1.25)]

cellModuleMatrix_norm <- as.data.frame(assay(pop, "logcounts"))

########## Subset factorized proportions matrix to just the modules of interest
popNorm2 <- cellModuleMatrix_norm[which(rownames(cellModuleMatrix_norm) %in% modules),]

########## Heatmap of individual celda module scores
popNorm3 <- t(scale(t(popNorm2)))

########## Reorder rows and columns & rename modules
popNorm3 <- popNorm3[c("L2", "L9", "L10", "L13", "L56", 
                       "L12", "L14", "L15", "L17", "L55",  
                       "L22", "L29",
                       "L23", "L24", "L25", "L32", "L59"),]

rownames(popNorm3)[which(rownames(popNorm3) == "L2")] <- "L2 (CD74)"
rownames(popNorm3)[which(rownames(popNorm3) == "L9")] <- "L9 (SCGB1A1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L10")] <- "L10 (SCGB3A1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L13")] <- "L13 (SLPI)"
rownames(popNorm3)[which(rownames(popNorm3) == "L56")] <- "L56 (MUC5B)"

rownames(popNorm3)[which(rownames(popNorm3) == "L12")] <- "L12 (BPIFB1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L14")] <- "L14 (MSMB)"
rownames(popNorm3)[which(rownames(popNorm3) == "L15")] <- "L15 (MUC5AC)"
rownames(popNorm3)[which(rownames(popNorm3) == "L17")] <- "L17 (TSPAN8)"
rownames(popNorm3)[which(rownames(popNorm3) == "L55")] <- "L55 (TFF3, LYZ)"

rownames(popNorm3)[which(rownames(popNorm3) == "L22")] <- "L22 (CEACAM5)"
rownames(popNorm3)[which(rownames(popNorm3) == "L29")] <- "L29 (TXN, SERPINB1)"

rownames(popNorm3)[which(rownames(popNorm3) == "L23")] <- "L23 (ADH1C)"
rownames(popNorm3)[which(rownames(popNorm3) == "L24")] <- "L24 (PSCA, MUC1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L25")] <- "L25 (LYPD2, PI3)"
rownames(popNorm3)[which(rownames(popNorm3) == "L32")] <- "L32 (AQP3, CYP26A1)"
rownames(popNorm3)[which(rownames(popNorm3) == "L59")] <- "L59 (S100A9)"

########## Prep heatmap annotations
column_ha <- HeatmapAnnotation(Cell_Type = pop$epithelial_CellType,
                               Smoking_Status = pop$Smoking_Status,
                               Histology = pop$Worst_Histology_at_TimePoint_Label,
                               SampleType = pop$SampleType,
                               # Celda_Cluster = pop$epithelial_celda_cluster,
                               col = list(Smoking_Status = smoking_colors,
                                          SampleType = sampletype_cols,
                                          Histology = histology_colors,
                                          Cell_Type = celltype_colors),
                               # Celda_Cluster = celda_cluster_cols),
                               show_annotation_name = TRUE,
                               show_legend = TRUE,
                               annotation_name_rot = 0)

row_ha <- rowAnnotation(Direction = c(rep("Club", 5), rep("Goblet", 5),
                                      rep("Perigoblet", 2), rep("Nasal Secretory", 5)),
                        col = list(Direction = c("Club" = "forestgreen", "Goblet" = "#3dd2e3", 
                                                 "Perigoblet" = "#fa9028", "Nasal Secretory" = "#cab2d6")),
                        show_annotation_name = FALSE,
                        show_legend = FALSE)

popNorm3[which(popNorm3 > 2)] <- 2
popNorm3[which(popNorm3 < -2)] <- -2

########## Plot & save heatmap
pdf(paste0("./Figures/Fig5D_DEModules_SecretoryCellTypes_Heatmap_Final.pdf"), width = 8, height = 8)
print(Heatmap(popNorm3,
              top_annotation = column_ha,
              right_annotation = row_ha,
              show_column_names = FALSE,
              column_split = pop$epithelial_CellType,
              show_heatmap_legend = TRUE,
              #column_title = NULL,
              column_title_rot = 45,
              column_title_gp = gpar(fontsize = 10),
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()


##################################### Figure 5E: Nasal secretory clusters vs each other
# nas_sec_cellprobs <- cellprobs[which(cellprobs$Cluster %in% c(27, 30, 34:37)),]
# diff_modules_nas_sec <- epi_clust_de(expr = nas_sec_cellprobs, clusters = unique(nas_sec_cellprobs$Cluster))
# diff_modules_res_nas_sec <- Reduce(diff_modules_nas_sec, f = rbind)
# diff_modules_res_significant_nas_sec <- diff_modules_res_nas_sec[which(diff_modules_res_nas_sec$SignificantPositive),]
# 
# saveRDS(diff_modules_res_nas_sec, "./Objects_Results/Fig5E_DEModules_NasalSecretoryClusters_Final.rds")
# saveRDS(diff_modules_res_significant_nas_sec, "./Objects_Results/Fig5E_DEModules_NasalSecretoryClusters_Significant_Final.rds")

diff_modules_res_nas_sec <- readRDS("./Objects_Results/Fig5E_DEModules_NasalSecretoryClusters_Final.rds")
diff_modules_res_significant_nas_sec <- readRDS("./Objects_Results/Fig5E_DEModules_NasalSecretoryClusters_Significant_Final.rds")

################### Individual clusters: Nasal secretory Heatmap (modules specific to each nas sec cluster)

########## Modules upregulated in nasal secretory clusters
modules <- unique(diff_modules_res_significant_nas_sec$Module[which(diff_modules_res_significant_nas_sec$Beta > 1)])

########## The factorizeMatrix function generates a "factorized" matrix between features/feature modules and samples/cell clusters/cells.
factoMat <- factorizeMatrix(sce, useAssay = "decontXcounts")
cellModuleMatrix <- factoMat$counts$cell

########## Standard normalization of counts
pop <- SingleCellExperiment(assays = SimpleList("factorizeCounts" = cellModuleMatrix), colData = colData(sce))
pop <- logNormCounts(pop, assay.type = "factorizeCounts", name = "logcounts")

pop <- pop[,which(pop$epithelial_celda_cluster %in% c(30, 36, 37, 34, 35, 27))] # 27, 30, 34:37

cellModuleMatrix_norm <- as.data.frame(assay(pop, "logcounts"))

########## Subset factorized proportions matrix to just the modules of interest
popNorm2 <- cellModuleMatrix_norm[which(rownames(cellModuleMatrix_norm) %in% modules),]

########## Heatmap of individual celda module scores
popNorm3 <- t(scale(t(popNorm2)))

########## Reorder rows and columns & rename modules
popNorm3 <- popNorm3[c("L70", 
                       "L26", 
                       "L9", "L11", "L12", "L14", "L10", 
                       "L15", "L17", "L55"),]

########## Prep heatmap annotations
column_ha <- HeatmapAnnotation(Smoking_Status = pop$Smoking_Status,
                               Histology = pop$Worst_Histology_at_TimePoint_Label,
                               # SampleType = pop$SampleType,
                               # Cell_Type = pop$epithelial_CellType,
                               # Celda_Cluster = pop$epithelial_celda_cluster,
                               col = list(Smoking_Status = smoking_colors,
                                          SampleType = sampletype_cols,
                                          Histology = histology_colors,
                                          Cell_Type = celltype_colors),
                               # Celda_Cluster = celda_cluster_cols),
                               show_annotation_name = TRUE,
                               show_legend = TRUE,
                               annotation_name_rot = 0)

#### Create a binary matrix for row annotations
num_modules <- length(unique(rownames(popNorm3)))
num_clusters <- 6
annotation_matrix <- matrix(0, nrow = num_modules, ncol = num_clusters)

rownames(annotation_matrix) <- unique(rownames(popNorm3))
colnames(annotation_matrix) <- c(30, 36, 37, 34, 35, 27)

# Fill in the matrix based on enrichment
for (i in 1:num_modules) {
  # Get the current module name
  module_num <- rownames(annotation_matrix)[i]
  
  # Find rows corresponding to the current module
  module_rows <- diff_modules_res_significant_nas_sec[diff_modules_res_significant_nas_sec$Module == module_num, ]
  
  # Extract the clusters as numeric values
  enriched_clusters <- as.numeric(as.character(module_rows$Cluster))
  
  # Filter out NA values
  enriched_clusters <- enriched_clusters[!is.na(enriched_clusters)]
  
  # Update the annotation matrix if enriched clusters exist
  if (length(enriched_clusters) > 0) {
    annotation_matrix[i, which(colnames(annotation_matrix) %in% enriched_clusters)] <- 1
  }
}

# Convert to a data frame for better handling
annotation_df <- as.data.frame(annotation_matrix)
colnames(annotation_df) <- paste0("Cluster_", colnames(annotation_df))

row_ha <- rowAnnotation(
  K30 = annotation_df$Cluster_30,
  K36 = annotation_df$Cluster_36,
  K37 = annotation_df$Cluster_37,
  K34 = annotation_df$Cluster_34,
  K35 = annotation_df$Cluster_35,
  K27 = annotation_df$Cluster_27,
  col = list(
    K30 = c("1" = "#dfc2f0", "0" = "white"),
    K36 = c("1" = "#8d36bf", "0" = "white"),
    K37 = c("1" = "#4a0870", "0" = "white"),
    K34 = c("1" = "#cd95ed", "0" = "white"),
    K35 = c("1" = "#ac5dd9", "0" = "white"),
    K27 = c("1" = "#e9ddf0", "0" = "white")
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "bottom",
  show_legend = FALSE # Suppress legend
)


popNorm3[which(popNorm3 > 2)] <- 2
popNorm3[which(popNorm3 < -2)] <- -2

pop$epithelial_celda_cluster <- unfactor(pop$epithelial_celda_cluster)
pop$epithelial_celda_cluster <- factor(pop$epithelial_celda_cluster, levels = c(30, 36, 37, 34, 35, 27))

########## Create & save heatmap
pdf(paste0("./Figures/Fig5E_DEModules_NasalSecretoryClusters_Heatmap_Final.pdf"), width = 9, height = 8)
print(Heatmap(popNorm3, 
              top_annotation = column_ha, 
              right_annotation = row_ha,
              show_column_names = FALSE, 
              column_split = pop$epithelial_celda_cluster,
              show_heatmap_legend = TRUE,
              #column_title = NULL,
              column_title_rot = 45,
              column_title_gp = gpar(fontsize = 10),
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()


##############################################################################################################
##################################### Figure 5F

sce <- sce_orig

########## Identify modules differentially expressed between clusters of the same cell type
f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")),
                                    colData = colData(sce))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")
sce_modular <- sce_modular[,-which(sce_modular$Smoking_Status == "Never")]
sce_modular <- sce_modular[,-which(sce_modular$SampleType == "Biopsy")]
sce_modular <- sce_modular[,-which(sce_modular$epithelial_CellType == "HG Basal Cells")]

cellprobs <- data.frame(t(assay(sce_modular,"LogNormalize")))
cellprobs$Cluster <- sce_modular$epithelial_celda_cluster ##
cellprobs$CellType <- sce_modular$epithelial_CellType
cellprobs$Sample <- sce_modular$Sample
cellprobs$Smoking_Status <- sce_modular$Smoking_Status
cellprobs$Smoking_Status <- factor(cellprobs$Smoking_Status, levels = c("Current", "Former"))
cellprobs$SampleType <- sce_modular$SampleType
cellprobs$SampleType <- factor(cellprobs$SampleType, levels = c("BronchialBrush", "NasalBrush"))
cellprobs$Site <- sce_modular$PCGA02_site

########## Model: module ~  sampletype + Smoking_Status + (1|sample)
epi_clust_de <- function(expr = cellprobs) {

    model_p <- data.frame(Module = names(expr)[1:70], 
                          P_Sm = rep(0, 70), 
                          Beta_Sm = rep(0, 70), 
                          # P_Br = rep(0, 70),
                          # Beta_Br = rep(0, 70),
                          P_Nas = rep(0, 70), 
                          Beta_Nas = rep(0, 70),
                          # P_SmBr = rep(0, 70),
                          # Beta_SmBr = rep(0, 70),
                          P_SmNas = rep(0, 70), 
                          Beta_SmNas = rep(0, 70),
                          # Group_Mean = rep(0, 70), 
                          # NotGroup_Mean = rep(0, 70), 
                          row.names = names(expr)[1:70]) # Smoking_P = rep(0, 55), Smoking_Beta = rep(0, 55), 
    
    for(i in rownames(model_p)) {
      
      model_construct<- as.formula(paste0(i," ~ Smoking_Status * SampleType"))
      model <- lme(model_construct,
                   random = ~ 1|Sample,
                   data = expr,
                   control = lmeControl(opt = "optim"),
                   # control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
                   na.action = na.omit)
      tTab <- summary(model)$tTable
      model_p[i,"Beta_Sm"] <- tTab["Smoking_StatusFormer","Value"]
      model_p[i,"P_Sm"] <- tTab["Smoking_StatusFormer","p-value"]
      model_p[i,"Beta_Nas"] <- tTab["SampleTypeNasalBrush","Value"]
      model_p[i,"P_Nas"] <- tTab["SampleTypeNasalBrush","p-value"]
      model_p[i,"Beta_SmNas"] <- tTab["Smoking_StatusFormer:SampleTypeNasalBrush","Value"]
      model_p[i,"P_SmNas"] <- tTab["Smoking_StatusFormer:SampleTypeNasalBrush","p-value"]
    }
    model_p$Q_Sm <- p.adjust(model_p$P_Sm, n = nrow(model_p), method = "fdr")
    model_p$Q_Nas <- p.adjust(model_p$P_Nas, n = nrow(model_p), method = "fdr")
    model_p$Q_SmNas <- p.adjust(model_p$P_SmNas, n = nrow(model_p), method = "fdr")
    model_p$LogQ_Sm <- -log10(model_p$Q_Sm)
    model_p$Significant <- ((model_p$Q_Sm < 0.05) & (abs(model_p$Beta_Sm) > 0.5))
    model_p$SignificantPositive <- ((model_p$Q_Sm < 0.05) & (model_p$Beta_Sm > 0.5))
    model_p <- model_p %>% arrange(desc(SignificantPositive), Q_Sm)
    model_p
}

################### Find modules up/down in smokers
# diff_modules_res <- epi_clust_de(expr = cellprobs)
# diff_modules_res_significant <- diff_modules_res[which(diff_modules_res$Significant),]

# saveRDS(diff_modules_res, "./Objects_Results/Fig5F_DEModules_Smoking_Final.rds")
# saveRDS(diff_modules_res_significant, "./Objects_Results/Fig5F_DEModules_Smoking_Significant_Final.rds")

diff_modules_res <- readRDS("./Objects_Results/Fig5F_DEModules_Smoking_Final.rds")
diff_modules_res_significant <- readRDS("./Objects_Results/Fig5F_DEModules_Smoking_Significant_Final.rds")

# write.xlsx(diff_modules_res_significant, "./Objects_Results/Fig5F_significant_model_res_Final.xlsx")

### modules <- diff_modules_res_significant$Module[which(diff_modules_res_significant$Cluster == "Current")]
modules <- diff_modules_res_significant$Module

########## Create & save violin plots
pdf(paste0("./Figures/Figure5F_L43_bySampleTypeSmokingStatus_Final.pdf"), height = 4.5, width = 4)
print(ggplot(cellprobs, aes(x = SampleType, y = L43, fill = Smoking_Status)) + 
        geom_violin(scale = "width") +
        theme_classic() + labs(title = "L43", y = "Expression", x = "") +
        scale_fill_manual(values = c("Former" = "grey40", "Current" = "red")) +
        theme(plot.title = element_text(hjust = 0.5), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              legend.position = "none")
)
dev.off()

pdf(paste0("./Figures/Figure5F_L52_bySampleTypeSmokingStatus_Final.pdf"), height = 4.5, width = 4)
print(ggplot(cellprobs, aes(x = SampleType, y = L52, fill = Smoking_Status)) + 
        geom_violin(scale = "width") +
        theme_classic() + labs(title = "L52", y = "Expression", x = "") +
        scale_fill_manual(values = c("Former" = "grey40", "Current" = "red")) +
        theme(plot.title = element_text(hjust = 0.5), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              legend.position = "none")
)
dev.off()


cellprobs$Smoking_Status <- factor(cellprobs$Smoking_Status, levels = c("Current", "Former"))
cellprobs_br <- cellprobs[which(cellprobs$SampleType == "BronchialBrush"),]
cellprobs_nas <- cellprobs[which(cellprobs$SampleType == "NasalBrush"),]

lmRes_interaction_L43 <- lme(L43 ~ Smoking_Status * SampleType, random = ~ 1|Sample, data = cellprobs, na.action = na.omit)
summary(lmRes_interaction_L43)
lmRes_interaction_L52 <- lme(L52 ~ Smoking_Status * SampleType, random = ~ 1|Sample, data = cellprobs, na.action = na.omit)
summary(lmRes_interaction_L52)

## Compare former smokers across sample types
emm_L43 <- emmeans(lmRes_interaction_L43, ~ SampleType | Smoking_Status)
pairs(emm_L43, by = "Smoking_Status")
emm_L52 <- emmeans(lmRes_interaction_L52, ~ SampleType | Smoking_Status)
pairs(emm_L52, by = "Smoking_Status")

##############################################################################################################
##################################### Figure 5G: DECAMP
##################### Read in data
t <- read.table("./072324_Biopsies+Brushes_noIonocytes_decontXcounts_zCellType_proteinCodingGenes_module_features_K40_L70.csv", sep = "\t")
t2 <- as.list(t)
t3 <- lapply(t2, function(z){ z[!is.na(z) & z != ""]})
names(t3) <- gsub("V", "L", names(t3))
modules <- diff_modules_res_significant$Module
smoking_modules <- t3[which(names(t3) %in% modules)]

#### DECAMP2
decamp2_annotation <- readRDS("/restricted/projectnb/decamp/annotation/20210805/20210805_DECAMP2_annotation.rds")


#### Bronch
decamp_bronch <- readRDS("/rprojectnb2/decamp/Analysis/COVID19_Analysis/GEO_Submission/geneExpression_bronch_raw_counts_20210318.rds")
decamp_bronch_annot <- readRDS("/rprojectnb2/decamp/Analysis/COVID19_Analysis/Bulk/DECAMP_BronchialBrushing/Data/DECAMP_Bronch_withCancerStats.rds")
cd_bronch <- as.data.frame(colData(decamp_bronch_annot))
rownames(cd_bronch) <- cd_bronch$kitnumber
cd_bronch <- cd_bronch[colnames(decamp_bronch),]

## Remove samples with no TIN median
no_TIN_median <- which(is.na(cd_bronch$TIN_median))
decamp_bronch <- decamp_bronch[,-no_TIN_median]
cd_bronch <- cd_bronch[-no_TIN_median,]

## Subset to DECAMP2
decamp_bronch <- decamp_bronch[,which(colnames(decamp_bronch) %in% decamp2_annotation$kitnumber)]
cd_bronch <- cd_bronch[which(cd_bronch$kitnumber %in% decamp2_annotation$kitnumber),]

#### Nasal
decamp_nasal <- readRDS("/rprojectnb2/decamp/Analysis/COVID19_Analysis/GEO_Submission/geneExpression_nasal_raw_counts_20210318.rds")
decamp_nasal_annot <- readRDS("/rprojectnb2/decamp/Analysis/COVID19_Analysis/Bulk/DECAMP_NasalBrushing/Data/combated_DECAMP_nasal_3batches_N288_050520.rds")
cd_nasal <- as.data.frame(colData(decamp_nasal_annot))
rownames(cd_nasal) <- cd_nasal$kitnumber
cd_nasal <- cd_nasal[colnames(decamp_nasal),]

## Subset to DECAMP2
decamp_nasal <- decamp_nasal[,which(colnames(decamp_nasal) %in% decamp2_annotation$kitnumber)]
cd_nasal <- cd_nasal[which(cd_nasal$kitnumber %in% decamp2_annotation$kitnumber),]

# #### Calculate residuals for decamp (already completed + saved)
# compute_residuals_limma <- function(count.matrix, tin, pheno, var1=NULL){
#   dge_BX <- DGEList(counts=count.matrix, samples = pheno)
#   keep.exprs <- filterByExpr(dge_BX)
#   # keep.exprs <- filterByExpr(dge_BX, group=dge_BX$samples[,var1])
#   # keep.exprs <- rowSums((edgeR::cpm(dge_BX$counts)) > 1) > ncol(dge_BX) *0.1
#   dge_BX <- dge_BX[keep.exprs,, keep.lib.sizes=FALSE]
#   dge_BX <- calcNormFactors(dge_BX, method = "TMM")
# 
#   design<-model.matrix(~as.numeric(tin))
# 
#   v<-voom(dge_BX,design=design)
#   fit<-lmFit(v,design)
#   fit<-eBayes(fit)
#   resid<-c()
#   resid<-residuals(fit,v)
#   return(resid)
# }
# 
# decamp_bronch.resid <- compute_residuals_limma(count.matrix = decamp_bronch,
#                                         tin = cd_bronch$TIN_median,
#                                         pheno = cd_bronch)
# 
# decamp_nasal.resid <- compute_residuals_limma(count.matrix = decamp_nasal,
#                                                tin = cd_nasal$TIN_median,
#                                                pheno = cd_nasal)

# ## Combine objects; subset to same gene list and then combine
# genes <- intersect(rownames(decamp_nasal.resid), rownames(decamp_bronch.resid))
# decamp_nasal.resid <- decamp_nasal.resid[which(rownames(decamp_nasal.resid) %in% genes),]
# decamp_bronch.resid <- decamp_bronch.resid[which(rownames(decamp_bronch.resid) %in% genes),]
# 
# saveRDS(decamp_nasal.resid, "./Objects_Results/Fig5G_DECAMP2_nasal_resid.rds")
# saveRDS(decamp_bronch.resid, "./Objects_Results/Fig5G_DECAMP2_bronch_resid.rds")
# 
# decamp.resid <- cbind(decamp_nasal.resid, decamp_bronch.resid)
# 
# saveRDS(decamp.resid, "./Objects_Results/Fig5G_DECAMP2_resid.rds")
# 
# Convert row names using R 3.6.0
# ensembl=useMart(host="www.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
# 
# ens <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),filters=c("ensembl_gene_id"),mart=ensembl,values=rownames(decamp.resid))
# decamp.resid2 <- decamp.resid[ens$ensembl_gene_id,]
# rownames(decamp.resid2) <- ens$hgnc_symbol
# 
# saveRDS(decamp.resid2, "./Objects_Results/Fig5G_DECAMP2_resid_symbols.rds")
# 
# # Read back in object after using R 3.6.0 for symbol conversion
# decamp.resid <- readRDS("./Objects_Results/Fig5G_DECAMP2_resid_symbols.rds")
# 
# ## Run GSVA
# decamp.gsva.scores <- gsva(expr = decamp.resid, gset.idx.list = smoking_modules, method = "gsva", kcdf = "Gaussian", mx.diff=1, abs.ranking=FALSE, verbose=TRUE)
# 
# ## save RDS file
# saveRDS(decamp.gsva.scores, "./Objects_Results/Fig5G_DECAMP2_GSVA_Scores.rds")

### Read in GSVA object
decamp.gsva.scores <- readRDS("./Objects_Results/Fig5G_DECAMP2_GSVA_Scores.rds")

## Add annotations
decamp.gsva.scores2 <- as.data.frame(t(decamp.gsva.scores))
decamp.gsva.scores2$TISSUE <- c(rep("Nasal", nrow(cd_nasal)), rep("Bronch", nrow(cd_bronch)))

decamp.gsva.scores2$SMOKING <- c(cd_nasal$smoking, cd_bronch$smk)
decamp.gsva.scores2$SMOKING <- gsub("Current smoker", "current", decamp.gsva.scores2$SMOKING)
decamp.gsva.scores2$SMOKING <- gsub("Former smoker", "former", decamp.gsva.scores2$SMOKING)

decamp.gsva.scores2$KITNUMBER <- c(cd_nasal$kitnumber, unfactor(cd_bronch$kitnumber))

# decamp.gsva.scores2 <- decamp.gsva.scores2[-which(decamp.gsva.scores2$SMOKING == "M"),]
decamp.gsva.scores2 <- decamp.gsva.scores2[-which(is.na(decamp.gsva.scores2$SMOKING)),]

decamp.gsva.scores2$SMOKING <- as.factor(decamp.gsva.scores2$SMOKING)
decamp.gsva.scores2$SMOKING <- relevel(decamp.gsva.scores2$SMOKING, ref = "former")

##### Visualization
decamp.gsva.resid.melt = reshape2::melt(decamp.gsva.scores2, id.vars = c("SMOKING", "TISSUE", "KITNUMBER"))

decamp.gsva.scores2$SMOKING <- relevel(decamp.gsva.scores2$SMOKING, ref = "current")

pdf(paste0("./Figures/Figure5G_L52_DECAMP2_GSVA_Scores_bySampleTypeSmokingStatus.pdf"), height = 4.5, width = 4)
ggplot(decamp.gsva.scores2, aes(x = TISSUE, y = L52, fill = SMOKING)) +
  geom_violin(scale = "width") +
  labs(title = "L52", y = "GSVA Score") +
  theme_classic() + 
  scale_fill_manual(values = c("red", "grey40"), name = "") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
dev.off()

pdf(paste0("./Figures/Figure5G_L43_DECAMP2_GSVA_Scores_bySampleTypeSmokingStatus.pdf"), height = 4.5, width = 4)
ggplot(decamp.gsva.scores2, aes(x = TISSUE, y = L43, fill = SMOKING)) +
  geom_violin(scale = "width") +
  labs(title = "L43", y = "GSVA Score") +
  theme_classic() + 
  scale_fill_manual(values = c("red", "grey40"), name = "") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
dev.off()

############## Modeling

decamp <- decamp.gsva.scores2 |>
  dplyr::mutate(
    TISSUE  = factor(TISSUE,  levels = c("Bronch","Nasal")),
    SMOKING = factor(SMOKING, levels = c("current","former")),
    KITNUMBER = factor(KITNUMBER)
  )

##### Interaction effects
m_decamp_L43 <- lme(L43 ~ SMOKING*TISSUE, random = ~ 1|KITNUMBER, data = decamp, na.action = na.omit)
summary(m_decamp_L43)

# report beta/p value
tt <- summary(m_decamp_L43)$tTable
beta_int_nas <- tt["SMOKINGformer:TISSUENasal","Value"]
p_int_nas <- tt["SMOKINGformer:TISSUENasal","p-value"]

m_decamp_L52 <- lme(L52 ~ SMOKING*TISSUE, random = ~ 1|KITNUMBER, data = decamp, na.action = na.omit)
summary(m_decamp_L52)

# report beta/p value
tt <- summary(m_decamp_L52)$tTable
beta_int_nas <- tt["SMOKINGformer:TISSUENasal","Value"]
p_int_nas <- tt["SMOKINGformer:TISSUENasal","p-value"]

###### Tissue types alone
em  <- emmeans(m_decamp_L43, ~ SMOKING | TISSUE)
eff <- contrast(em, list(Current_vs_Former=c(+1,-1)))
summary(eff, infer = c(TRUE, TRUE))

em  <- emmeans(m_decamp_L52, ~ SMOKING | TISSUE)
eff <- contrast(em, list(Current_vs_Former=c(+1,-1)))
summary(eff, infer = c(TRUE, TRUE))

## Compare former smokers across sample types
emm_L43 <- emmeans(m_decamp_L43, ~ TISSUE | SMOKING)
pairs(emm_L43, by = "SMOKING")
emm_L52 <- emmeans(m_decamp_L52, ~ TISSUE | SMOKING)
pairs(emm_L52, by = "SMOKING")

## Curr vs Former, averaged over tissues
emm_overall_L43 <- emmeans(m_decamp_L43, ~ SMOKING)  
pairs(emm_overall_L43)
emm_overall_L52 <- emmeans(m_decamp_L52, ~ SMOKING)   
pairs(emm_overall_L52)

##############################################################################################################
##################################### Figure 5H: Ponder
##################### Read in data
t <- read.table("./072324_Biopsies+Brushes_noIonocytes_decontXcounts_zCellType_proteinCodingGenes_module_features_K40_L70.csv", sep = "\t")
t2 <- as.list(t)
t3 <- lapply(t2, function(z){ z[!is.na(z) & z != ""]})
names(t3) <- gsub("V", "L", names(t3))
modules <- diff_modules_res_significant$Module
smoking_modules <- t3[which(names(t3) %in% modules)]

### Ponder data
ponder.symbol.dge <- readRDS("/restricted/projectnb/decamp/minyi/EGAC50000000169/COMBINED_OUTPUT/ponder.ensembl.combined.dge.rds")
ponder.symbol.bronchial.dge = ponder.symbol.dge[,ponder.symbol.dge$samples$Tissue == "Bronchial",]
ponder.symbol.nasal.dge = ponder.symbol.dge[,ponder.symbol.dge$samples$Tissue == "Nasal",]
ponder.symbol.nasal.dge$samples$cohort = ifelse(grepl("HV", ponder.symbol.nasal.dge$samples$CancerType), "Healthy", "Clinic")

ponder.symbol.nasal.HV.dge = ponder.symbol.nasal.dge[,ponder.symbol.nasal.dge$samples$cohort == "Healthy"]
ponder.symbol.nasal.CL.dge = ponder.symbol.nasal.dge[,ponder.symbol.nasal.dge$samples$cohort == "Clinic"]

####### Both nasal cohorts
# Removed former smokers who quit less than 12 months ago
ponder.symbol.nasal.dge = ponder.symbol.nasal.dge[, !(ponder.symbol.nasal.dge$samples$Smoking2 == "ex" & ponder.symbol.nasal.dge$samples$stop_months <= 12)]
ponder.symbol.nasal.dge = ponder.symbol.nasal.dge[,!(ponder.symbol.nasal.dge$samples$Smoking2 %in% c('1_to_12months', 'ex.less1month')),]

ponder.symbol.nasal.resid.design <- model.matrix(~Gender+Experiment+Age+RSeQC_TIN_median, data = ponder.symbol.nasal.dge$samples)
ponder.symbol.nasal.resid.voom <- voom(ponder.symbol.nasal.dge, ponder.symbol.nasal.resid.design, plot = F)
fit <- lmFit(ponder.symbol.nasal.resid.voom, ponder.symbol.nasal.resid.design)
fit<-eBayes(fit)
ponder.symbol.nasal.resid = residuals(fit,ponder.symbol.nasal.resid.voom)

####### Bronchial cohort
# Removed former smokers who quit less than 12 months ago
ponder.symbol.bronchial.dge = ponder.symbol.bronchial.dge[, !(ponder.symbol.bronchial.dge$samples$Smoking2 == "ex" & ponder.symbol.bronchial.dge$samples$stop_months <= 12)]
ponder.symbol.bronchial.dge = ponder.symbol.bronchial.dge[,!(ponder.symbol.bronchial.dge$samples$Smoking2 %in% c('1_to_12months', 'ex.less1month')),]

ponder.symbol.bronchial.resid.design <- model.matrix(~Gender+Experiment+Age+RSeQC_TIN_median, data = ponder.symbol.bronchial.dge$samples)
ponder.symbol.bronchial.resid.voom <- voom(ponder.symbol.bronchial.dge, ponder.symbol.bronchial.resid.design, plot = F)
fit <- lmFit(ponder.symbol.bronchial.resid.voom, ponder.symbol.bronchial.resid.design)
fit<-eBayes(fit)
ponder.symbol.bronchial.resid = residuals(fit,ponder.symbol.bronchial.resid.voom)

# ####### Bronchial & Nasal cohort, calculate GSVA together
# ponder.symbol.resid <- cbind(ponder.symbol.nasal.resid, ponder.symbol.bronchial.resid)
# 
# ### Convert ensembl to symbol
# LTP2_1 = readRDS("/restricted/projectnb/pulmseq/LTP2/BATCH_1/Output/Expression/LTP2_Gene_Expression.rds")
# 
# rownames(ponder.symbol.resid) <- rowData(LTP2_1)[match(rownames(ponder.symbol.resid), rownames(rowData(LTP2_1)),), "external_gene_name"]
# 
# Ponder.gsva.resid = gsva(expr = ponder.symbol.resid,
#                          gset.idx.list = smoking_modules,
#                          kcdf = "Gaussian",
#                          mx.diff = TRUE,
#                          abs.ranking = FALSE)
# 
# saveRDS(Ponder.gsva.resid, "./Objects_Results/Fig5H_Ponder_GSVA_Scores.rds")

Ponder.gsva.resid <- readRDS("./Objects_Results/Fig5H_Ponder_GSVA_Scores.rds")

Ponder.gsva.resid = as.data.frame(t(Ponder.gsva.resid))
Ponder.gsva.resid$SMOKING_STATUS = ponder.symbol.dge$samples[match(rownames(Ponder.gsva.resid),
                                                                   rownames(ponder.symbol.dge$samples)), "Smoking"]
Ponder.gsva.resid$TISSUE = ponder.symbol.dge$samples[match(rownames(Ponder.gsva.resid),
                                                           rownames(ponder.symbol.dge$samples)), "Tissue"]
Ponder.gsva.resid$SAMPLE = ponder.symbol.dge$samples[match(rownames(Ponder.gsva.resid),
                                                           rownames(ponder.symbol.dge$samples)), "SAMPLE_ID"]

Ponder.gsva.resid$PATIENT = ponder.symbol.dge$samples[match(rownames(Ponder.gsva.resid),
                                                            rownames(ponder.symbol.dge$samples)), "INDIVIDUAL_ID"]

Ponder.gsva.resid$SMOKING_STATUS = ifelse(Ponder.gsva.resid$SMOKING_STATUS == "current", "Current Smoker",
                                          ifelse(Ponder.gsva.resid$SMOKING_STATUS == "ex", "Former Smoker", "Never Smoker"))

Ponder.gsva.resid.melt = reshape::melt(Ponder.gsva.resid, id.vars = c("SMOKING_STATUS", "TISSUE", "PATIENT"))

Ponder.gsva.resid.melt$SMOKING_STATUS = factor(Ponder.gsva.resid.melt$SMOKING_STATUS, levels =c("Never Smoker","Former Smoker", "Current Smoker"))

smoking_colors <- readRDS("../smoking_colors.rds")
names(smoking_colors) <- c("Former Smoker", "Current Smoker", "Never Smoker")

##### Visualization
pdf(paste0("./Figures/Fig5G_L52_Ponder_GSVA_Scores_bySampleTypeSmokingStatus.pdf"), height = 4.5, width = 4)
ggplot(Ponder.gsva.resid, aes(x = TISSUE, y = L52, fill = SMOKING_STATUS)) +
  geom_violin(scale = "width") +
  labs(title = "L52", y = "GSVA Score") +
  theme_classic() + 
  scale_fill_manual(values = smoking_colors) +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
dev.off()

pdf(paste0("./Figures/Fig5G_L43_Ponder_GSVA_SCores_bySampleTypeSmokingStatus.pdf"), height = 4.5, width = 4)
ggplot(Ponder.gsva.resid, aes(x = TISSUE, y = L43, fill = SMOKING_STATUS)) +
  geom_violin(scale = "width") +
  labs(title = "L43", y = "GSVA Score") +
  theme_classic() + 
  scale_fill_manual(values = smoking_colors) +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
dev.off()

############## Modeling

ponder <- Ponder.gsva.resid |>
  dplyr::mutate(
    TISSUE = factor(TISSUE, levels = c("Bronchial","Nasal")),
    SMOKING_STATUS  = factor(SMOKING_STATUS, levels = c("Current Smoker","Former Smoker","Never Smoker")),
    PATIENT = factor(PATIENT)
  )

##### Interaction effects
m_ponder_L43 <- lme(L43 ~ SMOKING_STATUS*TISSUE, random = ~ 1|PATIENT, data = ponder, na.action = na.omit)
summary(m_ponder_L43)

tt <- summary(m_ponder_L43)$tTable
beta_int_FC <- tt["SMOKING_STATUSFormer Smoker:TISSUENasal","Value"]
p_int_FC <- tt["SMOKING_STATUSFormer Smoker:TISSUENasal","p-value"]
beta_int_NC <- tt["SMOKING_STATUSNever Smoker:TISSUENasal","Value"]
p_int_NC <- tt["SMOKING_STATUSNever Smoker:TISSUENasal","p-value"]

m_ponder_L52 <- lme(L52 ~ SMOKING_STATUS*TISSUE, random = ~ 1|PATIENT, data = ponder, na.action = na.omit)
summary(m_ponder_L52)

tt <- summary(m_ponder_L52)$tTable
beta_int_FC <- tt["SMOKING_STATUSFormer Smoker:TISSUENasal","Value"]
p_int_FC <- tt["SMOKING_STATUSFormer Smoker:TISSUENasal","p-value"]
beta_int_NC <- tt["SMOKING_STATUSNever Smoker:TISSUENasal","Value"]
p_int_NC <- tt["SMOKING_STATUSNever Smoker:TISSUENasal","p-value"]

###### Tissue types alone

# Current vs Former by tissue
em_L43 <- emmeans(m_ponder_L43, ~ SMOKING_STATUS | TISSUE)  # means by smoking within tissue
# Levels are c("Current Smoker","Former Smoker","Never Smoker")
cvf_L43 <- contrast(em_L43, list(Current_vs_Former = c(+1, -1, 0)))
summary(cvf_L43, infer = c(TRUE, TRUE))

em_L52 <- emmeans(m_ponder_L52, ~ SMOKING_STATUS | TISSUE)  # means by smoking within tissue
# Levels are c("Current Smoker","Former Smoker","Never Smoker")
cvf_L52 <- contrast(em_L52, list(Current_vs_Former = c(+1, -1, 0)))
summary(cvf_L52, infer = c(TRUE, TRUE))

# Current vs Never by tissue
em_L43 <- emmeans(m_ponder_L43, ~ SMOKING_STATUS | TISSUE)  # means by smoking within tissue
# Levels are c("Current Smoker","Former Smoker","Never Smoker")
cvn_L43 <- contrast(em_L43, list(Current_vs_Never = c(+1, 0, -1)))
summary(cvn_L43, infer = c(TRUE, TRUE))

em_L52 <- emmeans(m_ponder_L52, ~ SMOKING_STATUS | TISSUE)  # means by smoking within tissue
# Levels are c("Current Smoker","Former Smoker","Never Smoker")
cvn_L52 <- contrast(em_L52, list(Current_vs_Never = c(+1, 0, -1)))
summary(cvn_L52, infer = c(TRUE, TRUE))

###### Compare former smokers across sample types
emm_L43 <- emmeans(m_ponder_L43, ~ TISSUE | SMOKING_STATUS)
pairs(emm_L43, by = "SMOKING_STATUS")
emm_L52 <- emmeans(m_ponder_L52, ~ TISSUE | SMOKING_STATUS)
pairs(emm_L52, by = "SMOKING_STATUS")


## Curr vs Former, averaged over tissues
emm_overall_L43 <- emmeans(m_ponder_L43, ~ SMOKING_STATUS)  
pairs(emm_overall_L43)
emm_overall_L52 <- emmeans(m_ponder_L52, ~ SMOKING_STATUS)   
pairs(emm_overall_L52)
