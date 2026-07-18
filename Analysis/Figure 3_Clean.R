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
library(magrittr)
library(readxl)
library(jaccard)
library(GGally)
library(reshape2)
library(ggcorrplot)
library(glmmTMB)

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies")

sce <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/NonimmuneBiopsyCell_SCE_Final.rds")
basalgroup_cols <- c("Low Grade Basal Cells" = "blue", "High Grade Basal Cells" = "red")
hist_cols <- readRDS("histology_colors.rds")
smoke_cols <- readRDS("smoking_colors.rds")
cluster_colors <- distinctColors(n = nlevels(celdaClusters(sce)), hues = c("red", "orange", "yellow", "green", "blue", "purple", "brown"))
names(cluster_colors) <- levels(sce$KCluster)
basalonly_colors <- c(cluster_colors[1:11], rep("gray80",24))
names(basalonly_colors) <- levels(sce$KCluster)
hgb_basalcluster_colors <- basalonly_colors[c(paste0("K",1:6))]

f <- factorizeMatrix(sce, useAssay = "decontXcounts", type = "counts")
sce_modular <- SingleCellExperiment(assays = SimpleList("module_decontXcounts" = f$counts$cell),
                                    reducedDims = list("celda_UMAP" = reducedDim(altExp(sce), "celda_UMAP")), 
                                    colData = colData(sce))
sce_modular <- runSeuratNormalizeData(sce_modular, useAssay = "module_decontXcounts", normAssayName = "LogNormalize")

# Figure 3A: UMAP of Basal Cell Clusters
altExp(sce)$KCluster <- sce$KCluster
umap_coords <- reducedDim(altExp(sce), "celda_UMAP")
colnames(umap_coords) <- c("celda_UMAP_1", "celda_UMAP_2")
cluster_df <- cbind(data.frame(Clusters = sce$KCluster), umap_coords)
basalcluster_umap <- plotSCEDimReduceColData(inSCE = altExp(sce), reducedDimName = "celda_UMAP", conditionClass = "factor", colorBy = "KCluster", defaultTheme = FALSE, labelClusters = FALSE)
add_labels <- function(df, label_column, specificClusters, xcoord = "UMAP1", ycoord = "UMAP2", plt) {
  centroidList <- lapply(specificClusters, function(x) {
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
basalcluster_umap <- add_labels(cluster_df, "Clusters", specificClusters = paste0("K",1:11), xcoord = "celda_UMAP_1", ycoord = "celda_UMAP_2", plt = basalcluster_umap)
basalcluster_umap <- basalcluster_umap + theme_classic() + xlab("UMAP1") + ylab("UMAP2") + scale_color_manual(values = basalonly_colors) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "none")
pdf("BasalCluster_UMAP.pdf", width = 3, height = 3)
basalcluster_umap
dev.off()

## Compare module expression betwen HGB clusters
cellprobs <- data.frame(t(assay(sce_modular,"LogNormalize")))
cellprobs$Cluster <- sce$KCluster
mean_cellprobs <- cellprobs %>% group_by(Cluster) %>% summarize_all(mean)
mean_cellprobs <- data.frame(t(mean_cellprobs))
colnames(mean_cellprobs) <- mean_cellprobs["Cluster",]
mean_cellprobs <- mean_cellprobs[2:74,] 
mean_cellprobs <- sapply(mean_cellprobs, as.numeric)
mean_cellprobs <- data.frame(mean_cellprobs, row.names = colnames(cellprobs)[1:73])
cellprobs$CellType <- sce$CellType
cellprobs$Sample <- sce$Sample
cellprobs$Smoking_Status <- sce$Smoking_Status
cellprobs$Histology <- sce$Histology
cellprobs$PCGA02_Site <- sce$PCGA02_site

cluster_de <- function(expr = cellprobs, group_var = "Cluster", clusters) {
  expr_subset <- expr[which(expr[,group_var] %in% clusters),]
  each_model_res <- lapply(clusters, function(x) {
    expr_subset[,"GroupCluster"] <- ifelse(expr_subset[,group_var] == x, x, paste0("Not",x))
    expr_subset[,"GroupCluster"] <- factor(expr_subset[,"GroupCluster"],
                                           levels = c(paste0("Not",x),x),
                                           ordered = FALSE)
    expr_subset <- within(expr_subset, GroupCluster <- relevel(GroupCluster, ref = paste0("Not",x)))
    expr_subset[,"Smoking_Status"] <- factor(expr_subset[,"Smoking_Status"],
                                           levels = c("Former", "Current"),
                                           ordered = FALSE)
    expr_subset <- within(expr_subset, Smoking_Status <- relevel(Smoking_Status, ref = "Former"))
    expr_subset[,"PCGA02_Site"] <- factor(expr_subset[,"PCGA02_Site"],
                                             levels = c("Roswell", "UCL"),
                                             ordered = FALSE)
    expr_subset <- within(expr_subset, PCGA02_Site <- relevel(PCGA02_Site, ref = "Roswell"))
    model_p <- data.frame(Cluster = x, Module = names(expr_subset)[1:73], Group_Mean = rep(0, 73), NotGroup_Mean = rep(0, 73), Cluster_Beta = rep(0, 73), Cluster_P = rep(0, 73), Smoking_Beta = rep(0, 73), Smoking_P = rep(0, 73), row.names = names(expr_subset)[1:73])
    for(i in rownames(model_p)) {
      model_p[i,"Group_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] == x),i], na.rm = TRUE)
      model_p[i,"NotGroup_Mean"] <- mean(expr_subset[which(expr_subset[,group_var] != x),i], na.rm = TRUE)
      if(group_var == "Smoking_Status") {
        model_construct <- as.formula(paste0(i," ~ GroupCluster + PCGA02_Site"))
      } else {
        model_construct <- as.formula(paste0(i," ~ GroupCluster + Smoking_Status + PCGA02_Site"))
      }
      model <- lme(model_construct,
                   random = ~ 1|Sample,
                   data = expr_subset,
                   na.action = na.omit)
      tTab <- summary(model)$tTable
      model_var <- paste0("GroupCluster",x)
      model_p[i,"Cluster_Beta"] <- tTab[model_var,"Value"]
      model_p[i,"Cluster_P"] <- tTab[model_var,"p-value"]
      if(group_var != "Smoking_Status") {
        model_p[i,"Smoking_Beta"] <- tTab["Smoking_StatusCurrent","Value"]
        model_p[i,"Smoking_P"] <- tTab["Smoking_StatusCurrent","p-value"]
      }
    }
    model_p$Cluster_Q <- p.adjust(model_p$Cluster_P, n = nrow(model_p))
    model_p$Smoking_Q <- p.adjust(model_p$Smoking_P, n = nrow(model_p))
    model_p
  })
  names(each_model_res) <- clusters
  each_model_res
}
hgb_model_res <- cluster_de(clusters = paste0("K",1:6))
hgb_model_res <- Reduce(hgb_model_res,f = rbind)
hgb_model_res_sig <- hgb_model_res %>% filter(Cluster_Q < 0.05) %>% filter(Cluster_Beta > 0) %>% arrange(Cluster, Cluster_Q)
hgb_sigmods <- unique(hgb_model_res_sig$Module)
hgb_sigmods_ordered <- paste0("L",sort(as.numeric(gsub("L","",hgb_sigmods))))

## Figure 3B: HM of differentially expressed HGB modules
sce_basal <- sce_modular[, which(celdaClusters(sce) %in% c(1:6))]
cellAnnot <- data.frame("Cluster" = as.character(sce_basal$KCluster), "Sample" = sce_basal$Sample, "Smoking Status" = sce_basal$Smoking_Status, "Histology" = sce_basal$Histology, row.names = colnames(sce_basal), check.names = FALSE)
sample_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(sce_basal$Sample))); names(sample_cols) <- unique(cellAnnot$Sample)
CellAnnotColor <- list("Cluster" = hgb_basalcluster_colors, "Sample" = sample_cols, "Smoking Status" = smoke_cols, "Histology" = hist_cols) # see Figure 1.R for cluster_colors2, also it's saved

plotHM <- function(inSCE, assay = "LogNormalize", clusters, group_var, plot_var, features, cellAnnot, CellAnnotColor, show_legend) {
  cell_index <- which(inSCE[[group_var]] %in% clusters)
  sce_sub <- subsetSCECols(inSCE = inSCE, index = cell_index)
  plot_var_vec <- NULL
  if(is.character(plot_var)) {
    plot_var_vec <- sce_sub[[plot_var]]
  }
  f <- t(assay(sce_sub,assay))
  cell_probs_scale_mat <- t(scale(f[,features]))
  cell_probs_scale_mat[cell_probs_scale_mat > 2] <- 2; cell_probs_scale_mat[cell_probs_scale_mat < -2] <- -2
  colorScheme <- circlize::colorRamp2(c(min(cell_probs_scale_mat),(max(cell_probs_scale_mat) + min(cell_probs_scale_mat))/2, max(cell_probs_scale_mat)),
                                      c("blue", "white","red"))
  cat("Creating Heatmap Annotation\n")
  ca <- ComplexHeatmap::HeatmapAnnotation(df = cellAnnot, 
                                          col = CellAnnotColor, 
                                          show_legend = show_legend)
  cat("Plotting Heatmap\n")
  cell_probs_hm <- ComplexHeatmap::Heatmap(matrix = cell_probs_scale_mat, col = colorScheme, row_title = "Module", column_title = "Cell",
                                           cluster_rows = TRUE, cluster_columns = TRUE, column_split = plot_var_vec, 
                                           show_column_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, top_annotation = ca,
                                           heatmap_legend_param = list(title = "Scaled\nCell\nProbability"))
  cell_probs_hm <- draw(cell_probs_hm)
  col_order_list <- column_order(cell_probs_hm); row_order_list <- row_order(cell_probs_hm)
  list("heatmap" = cell_probs_hm, "column_order_list" = col_order_list, "row_order_list" = row_order_list, "matrix" = cell_probs_scale_mat) 
}
hgb_hm <- plotHM(inSCE = sce_modular, clusters = paste0("K",1:6), group_var = "KCluster", plot_var = "KCluster", features = hgb_sigmods_ordered, cellAnnot = cellAnnot, CellAnnotColor = CellAnnotColor, show_legend = c(TRUE, FALSE, FALSE, TRUE, TRUE))

## Figure 3C: Select differentially expressed HGB module violins
hgb_sigmods_subset_interesting <- c("L11","L15","L16","L30","L31","L43")
hgb_sigmods_subset_violin_interesting <- lapply(X = hgb_sigmods_subset_interesting, FUN = function(mod, BasalCluster = "Cluster") {
  cellprobs_sig <- cellprobs_hgb_sig[,c(mod, "Cluster")]
  sigmod_violin <- ggplot(cellprobs_sig, aes_string(x = BasalCluster, y = mod, fill = BasalCluster)) + geom_violin(scale = "width") +
    theme_classic() + labs(title = mod, y = "Expression", x = "") + scale_fill_manual(values = basalcluster_colors) +
    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_markdown(hjust = 0.5), legend.position = "none")
  sigmod_violin
})
pdf("hgb_sigmods_subset_violin_interesting.pdf", width = 6, height = 4)
plot_grid(plotlist = hgb_sigmods_subset_violin_interesting, nrow = 2, align = "hv")
dev.off()

## Compare outside module expression between HGB clusters 
vam_scores <- readRDS("Outside_Modules_VAM_Scores_All_Biopsy_Cells.rds") # See Figure 2_Clean.R
vam_scores <- vam_scores[colnames(sce),]
vam_scores$KCluster <- sce$KCluster
vam_scores$Smoking_Status <- sce$Smoking_Status
vam_scores$PCGA02_site <- sce$PCGA02_site
vam_scores$Sample <- sce$Sample
vam_scores_hgb <- vam_scores[which(vam_scores$KCluster %in% paste0("K",1:6)),]
vam_scores_hgb$KCluster <- factor(vam_scores_hgb$KCluster, levels = paste0("K",1:6), ordered = FALSE)
vam_scores_hgb <- within(vam_scores_hgb, KCluster <- relevel(KCluster, ref = "K1"))
vam_scores_hgb$Smoking_Status <- factor(vam_scores_hgb$Smoking_Status, levels = c("Former", "Current"), ordered = FALSE)
vam_scores_hgb <- within(vam_scores_hgb, Smoking_Status <- relevel(Smoking_Status, ref = "Former"))
vam_scores_hgb$PCGA02_site <- factor(vam_scores_hgb$PCGA02_site, levels = c("Roswell", "UCL"), ordered = FALSE)
vam_scores_hgb <- within(vam_scores_hgb, PCGA02_site <- relevel(PCGA02_site, ref = "Roswell"))

vam_model_logistic <- function(scores) {
  model_p <- data.frame(Module = rep(names(full_module_list), times = 6),
                        Cluster = rep(paste0("K",1:6), each = 28),
                        Cluster_Mean = rep(0, 28*6),
                        NotCluster_Mean = rep(0, 28*6),
                        Score_Logit = rep(0, 28*6),
                        Score_P = rep(0, 28*6),
                        Smoke_Logit = rep(0, 28*6),
                        Smoke_P = rep(0, 28*6))
  rownames(model_p) <- paste(model_p$Cluster, model_p$Module, sep = "_")
  for(i in 1:nrow(model_p)) {
    mod <- model_p$Module[i]
    scores_clustergrouped <- scores
    scores_clustergrouped$Group <- ifelse(scores_clustergrouped$KCluster == model_p$Cluster[i], "Group", "NotGroup")
    scores_clustergrouped$Group <- factor(scores_clustergrouped$Group, levels = c("NotGroup", "Group"), ordered = FALSE)
    scores_clustergrouped <- within(scores_clustergrouped, Group <- relevel(Group, ref = "NotGroup"))
    model_p[i, "Cluster_Mean"] <- mean(scores_clustergrouped[which(scores_clustergrouped$Group == "Group"), mod])
    model_p[i, "NotCluster_Mean"] <- mean(scores_clustergrouped[which(scores_clustergrouped$Group == "NotGroup"), mod])
    model_terms <- paste0("Group ~ ", mod, " + Smoking_Status + PCGA02_site + (1|Sample)")
    model_terms <- as.formula(model_terms)
    lmres <- glmmTMB(model_terms,
                     data=scores_clustergrouped,
                     family=binomial)
    coef <- summary(lmres)$coefficients$cond
    model_p[i,"Score_Logit"] <- coef[mod,"Estimate"]
    model_p[i,"Score_P"] <- coef[mod,"Pr(>|z|)"]
    model_p[i,"Smoke_Logit"] <- coef["Smoking_StatusCurrent","Estimate"]
    model_p[i,"Smoke_P"] <- coef["Smoking_StatusCurrent","Pr(>|z|)"]
  }
  model_p$Score_Q <- p.adjust(model_p$Score_P, method = "fdr")
  model_p$Smoke_Q <- p.adjust(model_p$Smoke_P, method = "fdr")
  list(model_results = model_p)
}
glmm_module_res_hgb <- vam_model_logistic(scores = vam_scores_hgb)
hgb_model_results <- glmm_module_res_hgb$model_results %>% arrange(Cluster, Score_Q)
scores_basal <- glmm_module_res_lgbvshgb$module_scores

## Figure 3D: Violins for Select Outside Modules in HGB Clusters
sig_diff_interesting <- c("Janes_cin", "Janes_epistroma_progup", "Janes_epistroma_progdown")
sig_diff_violin <- lapply(sig_diff_interesting, function(x) {
  plt <- ggplot(vam_scores_hgb, aes_string(x = "KCluster", y = x, fill = "KCluster")) + 
    geom_violin(scale = "width") + xlab("Cluster") + ylab("VAM Score") + ggtitle(x) + 
    theme_classic() + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + scale_fill_manual(values = hgb_basalcluster_colors) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")
})
sig_diff_interesting_violin <- plot_grid(plotlist = sig_diff_violin, align = "hv", ncol = 1)
pdf("Outside_Modules_HGB_Clusters.pdf", width = 3, height = 6)
sig_diff_interesting_violin
dev.off()

## Establish relationship between samples and clusters
## Figure 3E: SBP of clusters by sample
cellprobs_hgb_sig <- cellprobs[which(cellprobs$Cluster %in% paste0("K",1:6)),c(hgb_sigmods_ordered,"Cluster","Sample","Smoking_Status","Histology")]
# Use samples with at least 5 HGB
hgb_samp_table <- data.frame(table(cellprobs_hgb_sig$Sample))
has5hgb <- as.character(hgb_samp_table$Var1[hgb_samp_table$Freq >= 5])
cellprobs_hgb_sig <- cellprobs_hgb_sig[which(cellprobs_hgb_sig$Sample %in% has5hgb),]
sample_cluster_dist_df <- data.frame(table(cellprobs_hgb_sig$Sample, as.character(cellprobs_hgb_sig$Cluster)))
names(sample_cluster_dist_df) <- c("Sample","Cluster","Count")
sample_cluster_dist_df$Cluster <- as.character(sample_cluster_dist_df$Cluster)
sample_cluster_dist_df$Cluster <- factor(sample_cluster_dist_df$Cluster,
                                         levels = paste0("K",1:6),
                                         ordered = TRUE)
sample_cluster_dist_df$Count[which(sample_cluster_dist_df$Count == 0)] <- NA

cluster_sample_sbp <- sample_cluster_dist_df %>%
  group_by(Sample) %>%
  mutate(percentage = Count / sum(Count, na.rm = TRUE) * 100) %>%
  arrange(desc(percentage))
sbp_ordering <- c()
for(i in unique(cluster_sample_sbp$Sample)) {
  first_sample_index <- which(cluster_sample_sbp$Sample == i)[1]
  sbp_ordering <- c(sbp_ordering, first_sample_index)
}
cluster_sample_sbp_inorder <- cluster_sample_sbp[sbp_ordering,] %>% arrange(Cluster)
sample_cluster_dist_df$Sample <- factor(sample_cluster_dist_df$Sample,
                                        levels = cluster_sample_sbp_inorder$Sample,
                                        ordered = TRUE)

subject_sample <- merge(cluster_sample_sbp_inorder, unique(data.frame(Sample = sce$Sample, Subject = sce$Subject, TimePoint = sce$TimePoint, Location = sce$Anatomic_Site)))
m <- match(cluster_sample_sbp_inorder$Sample, subject_sample$Sample)
subject_sample <- subject_sample[m,]
subject_sample$Y <- 1

# Ordered by Subject, TP
subject_sample_tporder <- subject_sample %>% arrange(Subject, Location, TimePoint)
subject_sample_tporder$TimePoint[which(subject_sample_tporder$Sample == "54")] <- "T6"
subject_tp_order <- function(levels) {
  sample_cluster_dist_df$Sample <- factor(sample_cluster_dist_df$Sample,
                                          levels = levels,
                                          ordered = TRUE)
  sample_cluster_dist_sbp <- ggplot(sample_cluster_dist_df, aes(fill=Cluster, y=Count, x=Sample)) + 
    geom_bar(position="fill", stat="identity") + theme_classic() +
    labs(x = "Sample", y = "Cluster Percentage", title = "Distribution of Cells by Cluster per Sample") +
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
    scale_fill_manual(values = hgb_basalcluster_colors) +
    scale_x_discrete(expand = expansion(add = c(0,0))) + theme(plot.title = element_text(hjust = 0.5),
                                                                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  hist_df <- unique(cellprobs[,c("Sample","Histology")])
  m <- match(cluster_sample_sbp_inorder$Sample, hist_df$Sample)
  m <- m[!is.na(m)]
  hist_df <- hist_df[m,]
  hist_df$Sample <- factor(hist_df$Sample,
                           levels = levels(sample_cluster_dist_df$Sample),
                           ordered = TRUE)
  hist_df$X <- 1
  hist_samp_ggbar <- ggplot(hist_df, aes(x = Sample, y = X, fill = Histology)) + geom_bar(position = "fill", stat = "identity") +
    theme_minimal() + 
    scale_fill_manual(values = hist_cols) + xlab("") + ylab("Histology") +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          #axis.text = element_blank(), 
          legend.position = "right")
  
  smoke_df <- unique(cellprobs[,c("Sample","Smoking_Status")])
  m <- match(cluster_sample_sbp_inorder$Sample, smoke_df$Sample)
  m <- m[!is.na(m)]
  smoke_df <- smoke_df[m,]
  smoke_df$Sample <- factor(smoke_df$Sample,
                            levels = levels(sample_cluster_dist_df$Sample),
                            ordered = TRUE)
  smoke_df$X <- 1
  smoke_samp_ggbar <- ggplot(smoke_df, aes(x = Sample, y = X, fill = Smoking_Status)) + geom_bar(position = "fill", stat = "identity") +
    theme_minimal() + 
    scale_fill_manual(values = smoke_cols) + xlab("") + ylab("Smoking Status") +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          #axis.text = element_blank(), 
          legend.position = "right")
  
  predom_cluster <- data.frame(Cluster = cluster_sample_sbp_inorder[,"Cluster"], Sample = cluster_sample_sbp_inorder[,"Sample"], X = 1)
  predom_cluster$Sample <- as.character(predom_cluster$Sample)
  predom_cluster$Sample <- factor(predom_cluster$Sample,
                                  levels = levels,
                                  ordered = TRUE)
  m <- match(cluster_sample_sbp_inorder$Sample, predom_cluster$Sample)
  m <- m[!is.na(m)]
  predom_cluster_ggbar <- predom_cluster[m,]
  predom_cluster_ggbar <- ggplot(predom_cluster, aes(x = Sample, y = X, fill = Cluster)) + geom_bar(position = "fill", stat = "identity") +
    theme_minimal() + 
    scale_fill_manual(values = hgb_basalcluster_colors) + xlab("") + ylab("Predominant Cluster") +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          #axis.text = element_blank(), 
          legend.position = "none")
  
  subject_sample_tporder$ST <- paste(subject_sample_tporder$Subject, subject_sample_tporder$Location, subject_sample_tporder$TimePoint, sep = "_")
  subject_sample_tporder$Sample <- factor(subject_sample_tporder$Sample,
                                          levels = levels,
                                          ordered = TRUE)
  subject_text_plot <- ggplot(subject_sample_tporder, aes(x = Sample, y = Y, label = ST)) + geom_text() +
    theme_minimal() + 
    xlab("") + ylab("Subject") +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          #axis.text = element_blank(), 
          legend.position = "none")
  
  p <- plot_grid(plotlist = list(subject_text_plot, predom_cluster_ggbar, smoke_samp_ggbar, hist_samp_ggbar, sample_cluster_dist_sbp), align = "v", ncol = 1, rel_heights = c(1,1,1,1,4))
}
tporder_sample_sbp <- subject_tp_order(levels = subject_sample_tporder$Sample)
pdf("Sample_HGBCluster_SBP_wPredominantCluster+Subject_TPOrder_CorrectedSample54.pdf", width = 25, height = 11)
tporder_sample_sbp
dev.off()

########## ------------- SUPPLEMENTARY FIGURES/TABLES, FIGURE 3 ------------- ##########

## Supplementary Table 12: HGB 1-vs-all Module Results
openxlsx::write.xlsx(hgb_model_res, "HGB_DiffExp_Analysis_040225.xlsx")

## Supplementary Table 13: Outside Modules in HGB Clusters Module Results
openxlsx::write.xlsx(hgb_model_results, "Outside_Modules_VAM_LogR_HGB1vAll_Results.xlsx", rowNames = FALSE)

saveRDS(sce, "NonimmuneBiopsyCell_SCE_Final.rds")
