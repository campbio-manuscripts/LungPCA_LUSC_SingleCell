library(singleCellTK)
library(celda)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(magrittr)
library(scruff)
library(ggvenn)
library(forcats)
library(cowplot)
library(multipanelfigure)
library(kableExtra)
library(grid)
library(knitr)
library(RColorBrewer)
library(plotly)
library(pheatmap)
library(stringi)

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/")

sce <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/032624_Biopsies_Batches1-22_PCGenes_NoMTPseudo_Filtered_UMI>500_MT<0.4_p50<0.75.rds")

# Set up parameters for Celda
altExpName <- "featureSubset"
useAssay <- "decontXcounts"
reducedDimName <- "celda_UMAP"
maxFeatures <- 5000
initialK <- 5
maxK <- 60
initialL <- 5
maxL <- 125
markers <- c("EPCAM", "EGFR", "KRT4", "KRT6A","KRT6B", "KRT8", "KRT13", "KRT14", "KRT15", "KRT5", "TP63", "PDPN", "MAGEA4", "MAGEA6", "BCAM", "CDH1", "CLDN1",
             "SPDEF", "AGR2", "SCGB1A1", "SCGB3A1", "MUC5B",  "MUC5AC", "FOXJ1", "CAPS", "TPPP3", "TUBA4B", "TUBB4B", 
             "CFTR", "FOXI1", "ASCL3", "CALCA", "CHGA", "ASCL1",
             "COL1A1", "PDGFRA", "PDGFRB", "SPARC", "COL4A1", "COL4A2", "A2M", "MMP2", "SNAI1", "SNAI2", "ENG", "PECAM1", "VWF", "CA4", "CCN1", "ACTA2",
             "PTPRC", "CD3D", "CD4", "CCR7", "TCF7", "FOXP3", "CD8A", "GZMA", "GZMB", "PRF1", "NKG7", "PDCD1", "CTLA4", "LAG3", "TIGIT", 
             "CD19", "MS4A1", "CD79A", "MZB1", "LYZ", "CD14", "FCGR3A", "CD68", "CD163", "TYROBP", "ITGAM", "FCER1A", 
             "HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1",
             "MS4A2", "CPA3", "TPSAB1", "CSF3R", "FCGR3B", "S100A8", "S100A9", "IFITM1", "IFITM2", "IFITM3",
             "CEACAM5", "CEACAM6", "CYP1A1", "CYP1B1", "AKR1B10", "AKR1C1", "AKR1C2", "GSTP1", "GSTM1", "GSTM2", "GSTM3", "ALDH3A1", "NQO1", 
             "MKI67", "TOP2A", "MCM2", "MCM3", "MCM7", "CDK1", "CDKN2A", "PCNA", "CCNA1")
displayName <- "Gene_name"
exactMatch <- TRUE
showSetup <- TRUE
showSession <- TRUE
pdf <- TRUE

# List of cell annotations to include in the report
cellAnnot <- c("reads", "sum", "detected", "mt_percent", "percent_top_50", "ERCCcorrelation", "decontX_contamination", "PackYears")
# The following are plotted as plotly plots in the RMarkdown, and so don't have to be included in cellAnnotLabel: 
# Subject, Sample, Histology, Smoking_Status, Batch,  Flow
cellAnnotLabel <- c("Sex", "Age_Baseline", "Anatomic_Site", "PCGA02_site", "number_of_cells", "Dissociation_Enzyme", "Storage_Buffer",
                    "doubletFinder_doublet_label_resolution_1.5", "scDblFinder_doublet_call","scds_bcds_call","scds_cxds_call","scds_hybrid_call", 
                    "decontX_clusters")

## Decontaminate Counts Across Clusters
sce <- runDecontX(sce)

## Cluster genes into Modules
sce.subset <- sce[rowSums(round(decontXcounts(sce))) >= 1,]
dim(sce.subset) # 18348  9405
sce.subset <- selectFeatures(sce.subset, useAssay = useAssay)
varFilter <-
  ifelse(nrow(altExp(sce.subset, altExpName)) > maxFeatures, TRUE, FALSE)
if (varFilter) {
  # The HVG were already calculated in the sce object during QC (see Batch1-22_CellFiltering.R)
  o <-
    head(
      order(
        rowData(sce.subset)$seurat_variableFeatures_vst_varianceStandardized,
        decreasing = TRUE
      ),
      n = maxFeatures
    )
  altExp(sce.subset, altExpName) <- 
    subsetSCERows(sce.subset,
                  index = o,
                  returnAsAltExp = FALSE)
}

sampleLabel <- sce.subset$PCGA02_CellID_Template

moduleSplit <-
  recursiveSplitModule(
    sce.subset,
    initialL = initialL,
    maxL = maxL,
    sampleLabel = sampleLabel,
    altExpName = altExpName,
    useAssay = useAssay
  )

# Choose Module Number: 100
moduleSplitSelect <-
  subsetCeldaList(moduleSplit, params = list(L = 100))
initial.modules <- celdaModules(moduleSplitSelect)

## Cluster Cells into Clusters
cellSplit <-
  recursiveSplitCell(
    sce.subset,
    initialK = initialK,
    maxK = maxK,
    yInit = initial.modules,
    sampleLabel = sampleLabel,
    altExpName = altExpName,
    useAssay = useAssay)

# Choose Cluster Number: 30
sce.intermediate <- subsetCeldaList(cellSplit, params = list(K = 30))
sce.intermediate <- celdaUmap(sce.intermediate, useAssay = useAssay)
assay(sce.intermediate, "LogNormalizeCounts") <- assay(sce.intermediate, "LogNormalize")
sce.intermediate <- runSeuratNormalizeData(sce.intermediate, useAssay = "decontXcounts", normAssayName = "LogNormalize")
assay(altExp(sce.intermediate), "LogNormalizeCounts") <- assay(altExp(sce.intermediate), "LogNormalize")
altExp(sce.intermediate) <- runSeuratNormalizeData(altExp(sce.intermediate), useAssay = "decontXcounts", normAssayName = "LogNormalize")

sce.reorder <- sce.intermediate 
sce.reorder <- recodeClusterZ(sce.reorder, from = c(1:30), 
                              to = c(16,17,21,29,26,27,3,1,2,30,25,18,4,7,6,5,28,24,23,19,20,22,15,14,13,12,8,9,10,11))

## Split Modules based on expression across clusters
sce.split <- sce.reorder
sce.split <- splitModule(sce.split, module = 2, useAssay = useAssay)
sce.split <- splitModule(sce.split, module = 3, useAssay = useAssay)
sce.split <- splitModule(sce.split, module = 30, useAssay = useAssay)
sce.split <- splitModule(sce.split, module = 55, useAssay = useAssay)

# Verify stability of clusters and modules
sce.edit <- celda_CG(
  x = sce.split,
  useAssay = useAssay,
  altExpName = altExpName,
  sampleLabel = sce.split$PCGA02_LongID,
  K = 30,
  L = 104,
  alpha = 1,
  beta = 1,
  delta = 1,
  gamma = 1,
  stopIter = 500,
  maxIter = 500,
  splitOnIter = 10,
  splitOnLast = TRUE,
  seed = 12345,
  nchains = 1,
  zInitialize = "predefined",
  yInitialize = "predefined",
  countChecksum = NULL,
  zInit = celdaClusters(sce.split),
  yInit = celdaModules(sce.split),
  logfile = NULL,
  verbose = TRUE
)

sce.edit.reorder <- recodeClusterZ(sce.edit, from = 1:30,
                                   to = c(16,17,21,29,26,27,3,1,2,30,25,18,4,7,5,6,24,28,23,19,20,22,15,13,14,12,8,9,10,11))

## Render Visualization Report
moduleFileName_L104 <- paste0(format(Sys.Date(),"%m%d%y"),"_celda_K30_L104.tsv")
output_file_L104 <- paste0(format(Sys.Date(),"%m%d%y"),"_celda_K30_L104.html")
sceName_L104 <- paste0(format(Sys.Date(),"%m%d%y"),"_PCGA_K30_L104.rds")
output_dir <- "."
rmarkdown::render("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/Celda_CG_PlotResults.Rmd",
                  params = list(sce = sce.edit.reorder, 
                                sceName = sceName_L104,
                                altExpName = altExpName,
                                useAssay = useAssay,
                                reducedDimName = reducedDimName,
                                displayName = displayName,
                                cellAnnot = cellAnnot,
                                cellAnnotLabel = cellAnnotLabel,
                                exactMatch = exactMatch,
                                moduleFileName = moduleFileName_L104,
                                features = markers,
                                showSetup = showSetup,
                                showSession = showSession,
                                pdf= pdf),
                  output_file = output_file_L104,
                  output_dir = output_dir,
                  intermediates_dir = output_dir,
                  knit_root_dir = output_dir)

## Assign cell types based on expression of modules containing cell type-specific markers
sce.edit.reorder$CellType <- ""
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(1:7))] <- "Basal Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(8:9))] <- "Goblet Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(10:11))] <- "Club Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(12))] <- "SMG Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(13:15))] <- "Ciliated Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(16:17))] <- "Fibroblasts"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(18))] <- "Endothelial Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(19:21))] <- "CD8+ T Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(22:23))] <- "CD4+ T Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(24))] <- "B Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(25))] <- "Plasma Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(26:27))] <- "Macrophages"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(28))] <- "Dendritic Cells"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(29))] <- "Neutrophils"
sce.edit.reorder$CellType[which(celdaClusters(sce.edit.reorder) %in% c(30))] <- "Mast Cells"
altExp(sce.edit.reorder)$CellType <- sce.edit.reorder$CellType
saveRDS(sce.edit.reorder, paste0(format(Sys.Date(),"%m%d%y"),"_AllBiopsyCell_SCE_Final.rds"))
