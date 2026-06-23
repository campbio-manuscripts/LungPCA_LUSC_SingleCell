library(openxlsx)
library(SingleCellExperiment)
library(celda)

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies")
sce <- readRDS("121224_Nonimmune_SCE_BasalGroup.rds")
sce_meta <- readRDS("121224_Nonimmune_SCE_BasalGroup_PostMetadataRec_020325.rds")
sce_full <- readRDS("121224_AllCells_SCE_BasalGroup.rds")
sce_full_meta <- readRDS("121224_AllCells_SCE_BasalGroup_PostMetadataRec_020325.rds")
metadata_xlsx <- read.xlsx(xlsxFile = "../Data_Management/251125_PCGA2.0_scRNAseq_Biopsy+Brush_Metadata_2025-02-06_EEK.xlsx")

### All Cell Object Metadata Reconciliation
## Identifying unique columns of each object
setdiff(names(colData(sce_full)), names(colData(sce_full_meta))) # "PCGA02_LongID" "PCGA02_LongID_Batch"
setdiff(names(colData(sce_full_meta)), names(colData(sce_full)))
# "PCGA02_BiospecimenID" "PCGA02_BiospecimenID_Batch" "PCGA02_ParticipantID" "PCGA02_ParentID" "SampleType" 

# Finding Metdata Discrepancies
allcell_list <- list()
for(i in names(colData(sce_full_meta))) {
  if(i %in% names(colData(sce_full))) {
    allcell_list[[i]] <- which(sce_full_meta[[i]] != sce_full[[i]])
  }
}
## Columns for Reconciliation
sapply(allcell_list, length)
# experiment x
# PCGA02_CellID_Template x
# Batch x
# Sample x
# Flow x
# TimePoint x
# Age_Baseline x
# Anatomic_Site x
# Cohort x
# Subject_TimePoint x
# PackYears x
# Ethnicity x
# Race x

View(unique(data.frame(sce_full_meta$Sample[allcell_list[["experiment"]]], Meta = sce_full_meta$experiment[allcell_list[["experiment"]]], Orig = sce_full$experiment[allcell_list[["experiment"]]])))
# Keeping Meta
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["PCGA02_CellID_Template"]]], Meta = sce_full_meta$PCGA02_CellID_Template[allcell_list[["PCGA02_CellID_Template"]]], Orig = sce_full$PCGA02_CellID_Template[allcell_list[["PCGA02_CellID_Template"]]])))
# Keeping Meta   
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["Batch"]]], Meta = sce_full_meta$Batch[allcell_list[["Batch"]]], Orig = sce_full$Batch[allcell_list[["Batch"]]])))
batch_orig <- as.character(sce_full$Batch)
sce_full_meta$Batch[allcell_list[["Batch"]]] <- batch_orig[allcell_list[["Batch"]]] # Switching from Orig
sce_full_meta$Batch <- factor(sce_full_meta$Batch, levels = levels(sce_full$Batch), ordered = TRUE)
View(unique(data.frame(Meta = sce_full_meta$Sample[allcell_list[["Sample"]]], Orig = sce_full$Sample[allcell_list[["Sample"]]])))
# Keeping Meta   
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["Flow"]]], Meta = sce_full_meta$Flow[allcell_list[["Flow"]]], Orig = sce_full$Flow[allcell_list[["Flow"]]])))
# Keeping Meta   
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["TimePoint"]]], Meta = sce_full_meta$TimePoint[allcell_list[["TimePoint"]]], Orig = sce_full$TimePoint[allcell_list[["TimePoint"]]])))
# Keeping Meta 
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["Age_Baseline"]]], Meta = sce_full_meta$Age_Baseline[allcell_list[["Age_Baseline"]]], Orig = sce_full$Age_Baseline[allcell_list[["Age_Baseline"]]])))
# Keeping Meta 
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["Anatomic_Site"]]], Meta = sce_full_meta$Anatomic_Site[allcell_list[["Anatomic_Site"]]], Orig = sce_full$Anatomic_Site[allcell_list[["Anatomic_Site"]]])))
# Keeping Meta 
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["Cohort"]]], Meta = sce_full_meta$Cohort[allcell_list[["Cohort"]]], Orig = sce_full$Cohort[allcell_list[["Cohort"]]])))
sce_full_meta$Cohort <- gsub("Cross-sectional", "Cross-Sectional", sce_full_meta$Cohort) # Keeping Meta 
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["Subject_TimePoint"]]], Meta = sce_full_meta$Subject_TimePoint[allcell_list[["Subject_TimePoint"]]], Orig = sce_full$Subject_TimePoint[allcell_list[["Subject_TimePoint"]]])))
# Keeping Meta 
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["PackYears"]]], Meta = sce_full_meta$PackYears[allcell_list[["PackYears"]]], Orig = sce_full$PackYears[allcell_list[["PackYears"]]])))
# Keeping Meta 
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["Ethnicity"]]], Meta = sce_full_meta$Ethnicity[allcell_list[["Ethnicity"]]], Orig = sce_full$Ethnicity[allcell_list[["Ethnicity"]]])))
# Keeping Meta 
View(unique(data.frame(sce_full_meta$Sample[allcell_list[["Race"]]], Meta = sce_full_meta$Race[allcell_list[["Race"]]], Orig = sce_full$Race[allcell_list[["Race"]]])))
# Keeping Meta 

# Copy correct metadata to altExp
for(i in names(colData(sce_full_meta))) {
  if(i %in% names(colData(altExp(sce_full_meta)))) {
    altExp(sce_full_meta)[[i]] <- sce_full_meta[[i]]
  }
}

rm(batch_orig, i)

### Nonimmune Cell Object Metadata Harmonization
## Identifying unique columns of each object
setdiff(names(colData(sce)), names(colData(sce_meta))) # "PCGA02_LongID" "HistologyCondense"
setdiff(names(colData(sce_meta)), names(colData(sce)))
# "PCGA02_BiospecimenID" "PCGA02_ParticipantID" "PCGA02_ParentID" "SampleType" "PCGA02_BiospecimenID_Batch"

noni_list <- list()
for(i in names(colData(sce_meta))) {
  if(i %in% names(colData(sce))) {
    noni_list[[i]] <- which(sce_meta[[i]] != sce[[i]])
  }
}
## Columns for Reconciliation
sapply(noni_list, length)
# experiment x
# PCGA02_CellID_Template x
# Sample x
# Flow x
# TimePoint x
# Age_Baseline x
# Anatomic_Site x
# Cohort x
# Subject_TimePoint x
# PackYears x
# Ethnicity x
# Race x
# CellType x

View(unique(data.frame(sce_meta$Sample[noni_list[["experiment"]]], Meta = sce_meta$experiment[noni_list[["experiment"]]], Orig = sce$experiment[noni_list[["experiment"]]])))
# Keeping Meta
View(unique(data.frame(sce_meta$Sample[noni_list[["PCGA02_CellID_Template"]]], Meta = sce_meta$PCGA02_CellID_Template[noni_list[["PCGA02_CellID_Template"]]], Orig = sce$PCGA02_CellID_Template[noni_list[["PCGA02_CellID_Template"]]])))
# Keeping Meta   
View(unique(data.frame(sce_meta$Sample[noni_list[["Batch"]]], Meta = sce_meta$Batch[noni_list[["Batch"]]], Orig = sce$Batch[noni_list[["Batch"]]])))
batch_orig <- as.character(sce$Batch)
sce_meta$Batch[noni_list[["Batch"]]] <- batch_orig[noni_list[["Batch"]]] # Switching from Orig
sce_meta$Batch <- factor(sce_meta$Batch, levels = levels(sce$Batch), ordered = TRUE)
View(unique(data.frame(Meta = sce_meta$Sample[noni_list[["Sample"]]], Orig = sce$Sample[noni_list[["Sample"]]])))
# Keeping Meta   
View(unique(data.frame(sce_meta$Sample[noni_list[["Flow"]]], Meta = sce_meta$Flow[noni_list[["Flow"]]], Orig = sce$Flow[noni_list[["Flow"]]])))
# Keeping Meta   
View(unique(data.frame(sce_meta$Sample[noni_list[["TimePoint"]]], Meta = sce_meta$TimePoint[noni_list[["TimePoint"]]], Orig = sce$TimePoint[noni_list[["TimePoint"]]])))
# Keeping Meta 
View(unique(data.frame(sce_meta$Sample[noni_list[["Age_Baseline"]]], Meta = sce_meta$Age_Baseline[noni_list[["Age_Baseline"]]], Orig = sce$Age_Baseline[noni_list[["Age_Baseline"]]])))
# Keeping Meta 
View(unique(data.frame(sce_meta$Sample[noni_list[["Anatomic_Site"]]], Meta = sce_meta$Anatomic_Site[noni_list[["Anatomic_Site"]]], Orig = sce$Anatomic_Site[noni_list[["Anatomic_Site"]]])))
# Keeping Meta 
View(unique(data.frame(sce_meta$Sample[noni_list[["Cohort"]]], Meta = sce_meta$Cohort[noni_list[["Cohort"]]], Orig = sce$Cohort[noni_list[["Cohort"]]])))
sce_meta$Cohort <- gsub("Cross-sectional", "Cross-Sectional", sce_meta$Cohort) # Keeping Meta 
View(unique(data.frame(sce_meta$Sample[noni_list[["Subject_TimePoint"]]], Meta = sce_meta$Subject_TimePoint[noni_list[["Subject_TimePoint"]]], Orig = sce$Subject_TimePoint[noni_list[["Subject_TimePoint"]]])))
# Keeping Meta 
View(unique(data.frame(sce_meta$Sample[noni_list[["PackYears"]]], Meta = sce_meta$PackYears[noni_list[["PackYears"]]], Orig = sce$PackYears[noni_list[["PackYears"]]])))
# Keeping Meta 
View(unique(data.frame(sce_meta$Sample[noni_list[["Ethnicity"]]], Meta = sce_meta$Ethnicity[noni_list[["Ethnicity"]]], Orig = sce$Ethnicity[noni_list[["Ethnicity"]]])))
# Keeping Meta 
View(unique(data.frame(sce_meta$Sample[noni_list[["Race"]]], Meta = sce_meta$Race[noni_list[["Race"]]], Orig = sce$Race[noni_list[["Race"]]])))
# Keeping Meta 
View(unique(data.frame(sce_meta$Sample[noni_list[["CellType"]]], Meta = sce_meta$CellType[noni_list[["CellType"]]], Orig = sce$CellType[noni_list[["CellType"]]])))
sce_meta$CellType <- sce$CellType

# Copy correct metadata to altExp
for(i in names(colData(sce_meta))) {
  if(i %in% names(colData(altExp(sce_meta)))) {
    altExp(sce_meta)[[i]] <- sce_meta[[i]]
  }
}

saveRDS(sce_full_meta, "AllBiopsyCell_SCE_Final.rds")
saveRDS(sce_meta, "NonimmuneBiopsyCell_SCE_Final.rds")
