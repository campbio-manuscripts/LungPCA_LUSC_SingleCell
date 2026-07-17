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
library(multipanelfigure)
library(ggtext)
library(forcats)
library(openxlsx)

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/")

sce_full <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/070224_AllBiopsyCell_SCE_Final.rds")

nonimmune <- which(celdaClusters(sce_full) %in% c(1:18))
sce <- subsetSCECols(sce_full, index = nonimmune)
altExp(sce) <- c()

# Set up parameters for Celda
altExpName <- "featureSubset"
useAssay <- "decontXcounts"
reducedDimName <- "celda_UMAP"
maxFeatures <- 5000
initialK <- 1
maxK <- 60
initialL <- 1
maxL <- 125
markers <- c("EPCAM", "EGFR", "KRT4", "KRT6A","KRT6B", "KRT8", "KRT13", "KRT14", "KRT15", "KRT5", "TP63", "PDPN", "MAGEA4", "MAGEA6", "BCAM", "CDH1", "CLDN1",
             "SPDEF", "AGR2", "SCGB1A1", "SCGB3A1", "MUC5B",  "MUC5AC", "FOXJ1", "CAPS", "TPPP3", "TUBA4B", "TUBB4B", "PRB1", "PRB2", "PRB3",
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

## Choose features that are variable across nonimmune cells, rather than all cells
sce.subset <- sce[rowSums(round(decontXcounts(sce))) >= 1,]
dim(sce.subset) # 17719  4579
sce.subset <- selectFeatures(sce.subset, useAssay = useAssay)
dim(altExp(sce.subset)) # 9560 4579
varFilter <-
  ifelse(nrow(altExp(sce.subset, altExpName)) > maxFeatures, TRUE, FALSE)
metadata(sce.subset)$seurat$obj <- c() # otherwise, the gene variance calculation from sce_full does not get overwritten and throws an error
if (varFilter) {
  temp.sce <- sce.subset
  seuratObject <- convertSCEToSeurat(temp.sce, countsAssay = useAssay)
  seuratObject <- Seurat::FindVariableFeatures(seuratObject, 
                                               selection.method = "vst", nfeatures = 2000, verbose = TRUE)
  temp.sce <- singleCellTK:::.addSeuratToMetaDataSCE(temp.sce, seuratObject)
  rowData(temp.sce)$seurat_variableFeatures_vst_varianceStandardized <- methods::slot(temp.sce@metadata$seurat$obj, 
                                                                                      "assays")[["RNA"]]@meta.features$vst.variance.standardized
  temp.sce <-
    runSeuratFindHVG(inSCE = temp.sce, useAssay = useAssay)
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

highexpr_genes <- which(rowSums(decontXcounts(altExp(sce.subset))) >= 1)
altExp(sce.subset) <- subsetSCERows(altExp(sce.subset), index = highexpr_genes, returnAsAltExp = FALSE)
sampleLabel <- sce.subset$PCGA02_CellID_Template

## Cluster genes into Modules
moduleSplit <-
  recursiveSplitModule(
    sce.subset,
    initialL = initialL,
    maxL = maxL,
    sampleLabel = sampleLabel,
    altExpName = altExpName,
    useAssay = useAssay
  )

# Choose Module Number: 70
moduleSplitSelect <-
  subsetCeldaList(moduleSplit, params = list(L = 70))
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

# Choose Cluster Number: 35
sce.intermediate <- subsetCeldaList(cellSplit, params = list(K = 35))
sce.intermediate <- celdaUmap(sce.intermediate, useAssay = useAssay)
sceName <- paste0(format(Sys.Date(),"%m%d%y"),"_Biopsies_Batches1-22_Nonimmune_Celda_K35_L70.rds")

## Split Modules based on expression across clusters
sce.split <- sce.intermediate
sce.split <- splitModule(sce.split, module = 58, n = 3, useAssay = useAssay)
sce.split <- splitModule(sce.split, module = 59, useAssay = useAssay)

# Verify stability of clusters and modules
sce.edit <- celda_CG(
  x = sce.split,
  useAssay = useAssay,
  altExpName = altExpName,
  sampleLabel = sce.split$PCGA02_LongID,
  K = 35,
  L = 73,
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

sce.edit <- celdaUmap(sce.edit, useAssay = useAssay)
sce.reorder <- sce.edit
sce.reorder <- recodeClusterZ(sce.reorder, from = c(1:35), 
                              to = c(24,23,19,18,8,12,7,11,9,10,17,20,21,22,16,13,15,14,34,35,30,28,27,25,26,29,31,33,32,6,3,1,4,5,2))

## Render Visualization Report
moduleFileName_L73_reorder <- paste0(format(Sys.Date(),"%m%d%y"),"_Nonimmune_celda_K35_L73_CeldaEdit.tsv")
output_file_L73_reorder <- paste0(format(Sys.Date(),"%m%d%y"),"_Nonimmune_celda_K35_L73_CeldaEdit.html")
sceName_L73_reorder <- paste0(format(Sys.Date(),"%m%d%y"),"_Nonimmune_celda_K35_L73_CeldaEdit_CellTypeAnnotated.rds")
output_dir <- "."
rmarkdown::render("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/Celda_CG_PlotResults.Rmd",
                  params = list(sce = sce.reorder, 
                                sceName = sceName_L73_reorder,
                                altExpName = altExpName,
                                useAssay = useAssay,
                                reducedDimName = reducedDimName,
                                displayName = displayName,
                                cellAnnot = cellAnnot,
                                cellAnnotLabel = cellAnnotLabel,
                                exactMatch = exactMatch,
                                moduleFileName = moduleFileName_L73_reorder,
                                features = nonimmune_markers,
                                showSetup = showSetup,
                                showSession = showSession,
                                pdf= pdf),
                  output_file = output_file_L73_reorder,
                  output_dir = output_dir,
                  intermediates_dir = output_dir,
                  knit_root_dir = output_dir)

## Reassign cell types based on expression of nonimmune modules 
sce.reorder$CellTypeAllCells <- sce.reorder$CellType
altExp(sce.reorder)$CellTypeAllCells <- sce.reorder$CellTypeAllCells
sce.reorder$CellType <- ""
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(1:9))] <- "Basal Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(10))] <- "KRT5+/SCGB1A1+ Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(11))] <- "KRT5+/MUC5B+ Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(12))] <- "Perigoblet Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(13:17))] <- "Goblet Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(18:22))] <- "Club Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(23))] <- "Mucous SMG Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(24))] <- "Serous SMG Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(25:30))] <- "Ciliated Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(31))] <- "Airway Smooth Muscle Cells"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(32:33))] <- "Fibroblasts"
sce.reorder$CellType[which(celdaClusters(sce.reorder) %in% c(34:35))] <- "Endothelial Cells"
altExp(sce.reorder)$CellType <- sce.reorder$CellType
saveRDS(sce.reorder, paste0(format(Sys.Date(),"%m%d%y"),"_NonimmuneBiopsyCell_SCE_Final.rds"))

