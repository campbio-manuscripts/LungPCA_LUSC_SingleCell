library(openxlsx)
library(SingleCellExperiment)
library(celda)

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies")

metadata_xlsx <- read.xlsx("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Data_Management/PCGA2.0_scRNAseq_Biopsy+Brush_Metadata_2024-12_EEK.xlsx")
sce_full <- readRDS("121224_AllCells_SCE_BasalGroup.rds")
sce_full$PCGA02_ParticipantID <- ""
coldata_sce_full <- data.frame(colData(sce_full))
sce <- readRDS("121224_Nonimmune_SCE_BasalGroup.rds")
sce$PCGA02_ParticipantID <- ""
coldata_sce <- data.frame(colData(sce))
### 02/19/25: This will only happen after some additional metadata reconciliation
saveRDS(sce_full, "121224_AllCells_SCE_BasalGroup_PreMetadataRec_020325.rds")
saveRDS(sce, "121224_Nonimmune_SCE_BasalGroup_PreMetadataRec_020325.rds")

## Subject Level Metadata
subject_metadata_xlsx <- metadata_xlsx[,c("PCGA02_ParticipantID", "Subject", "PCGA02_site", "Age_Baseline_DSB", "Sex", "Smoking_Status", "Cohort", "PackYears", "Ethnicity", "Race")]
subject_metadata_xlsx <- unique(subject_metadata_xlsx)
subject_metadata_xlsx$Smoking_Status <- gsub(" smoker","", subject_metadata_xlsx$Smoking_Status)
coldata_sce_full_subject <- coldata_sce_full[,c("PCGA02_ParticipantID", "Subject", "PCGA02_site", "Age_Baseline", "Sex", "Smoking_Status", "Cohort", "PackYears", "Ethnicity", "Race")]
coldata_sce_full_subject <- unique(coldata_sce_full_subject)
colnames(coldata_sce_full_subject)[4] <- "Age_Baseline_DSB"
coldata_sce_subject <- coldata_sce[,c("PCGA02_ParticipantID", "Subject", "PCGA02_site", "Age_Baseline", "Sex", "Smoking_Status", "Cohort", "PackYears", "Ethnicity", "Race")]
coldata_sce_subject <- unique(coldata_sce_subject)
colnames(coldata_sce_subject)[4] <- "Age_Baseline_DSB"

full_subject_reconciliation <- list()
for(i in subject_metadata_xlsx$Subject) {
  if(i %in% unique(coldata_sce_full_subject$Subject)) {
    xlsx_data <- subject_metadata_xlsx[which(subject_metadata_xlsx$Subject == i),]
    sce_full_data <- coldata_sce_full_subject[which(coldata_sce_full_subject$Subject == i),]
    full_subject_reconciliation[[i]] <- rbind(xlsx_data, sce_full_data)
  }
}
full_subject_reconciliation <- full_subject_reconciliation[!sapply(full_subject_reconciliation,is.null)]

# No discrepancies or missing values for PCGA02_site, Sex, or Smoking_Status
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["PCGA02_ParticipantID"]]))}))
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["PCGA02_site"]]))}))
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["Age_Baseline_DSB"]]))}))
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["Sex"]]))}))
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["Smoking_Status"]]))}))
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["Cohort"]]))}))
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["PackYears"]]))}))
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["Ethnicity"]]))}))
table(sapply(full_subject_reconciliation, function(x) {length(unique(x[["Race"]]))}))
for(i in subject_metadata_xlsx$Subject) {
  if(i %in% unique(sce_full$Subject)) {
    sce_m <- which(sce_full$Subject == i)
    meta_m <- which(subject_metadata_xlsx$Subject == i)
    sce_full$PCGA02_ParticipantID[sce_m] <- subject_metadata_xlsx$PCGA02_ParticipantID[meta_m]
    sce_full$Age_Baseline[sce_m] <- subject_metadata_xlsx$Age_Baseline_DSB[meta_m]
    sce_full$Cohort[sce_m] <- subject_metadata_xlsx$Cohort[meta_m]
    sce_full$PackYears[sce_m] <- subject_metadata_xlsx$PackYears[meta_m]
    sce_full$Ethnicity[sce_m] <- subject_metadata_xlsx$Ethnicity[meta_m]
    sce_full$Race[sce_m] <- subject_metadata_xlsx$Race[meta_m]
  }
}

nonimmune_subject_reconciliation <- list()
for(i in subject_metadata_xlsx$Subject) {
  if(i %in% unique(coldata_sce_subject$Subject)) {
    xlsx_data <- subject_metadata_xlsx[which(subject_metadata_xlsx$Subject == i),]
    sce_full_data <- coldata_sce_subject[which(coldata_sce_subject$Subject == i),]
    nonimmune_subject_reconciliation[[i]] <- rbind(xlsx_data, sce_full_data)
  }
}
nonimmune_subject_reconciliation <- nonimmune_subject_reconciliation[!sapply(nonimmune_subject_reconciliation,is.null)]

# No discrepancies or missing values for PCGA02_site, Sex, or Smoking_Status
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["PCGA02_ParticipantID"]]))}))
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["PCGA02_site"]]))}))
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["Age_Baseline_DSB"]]))}))
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["Sex"]]))}))
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["Smoking_Status"]]))}))
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["Cohort"]]))}))
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["PackYears"]]))}))
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["Ethnicity"]]))}))
table(sapply(nonimmune_subject_reconciliation, function(x) {length(unique(x[["Race"]]))}))
for(i in subject_metadata_xlsx$Subject) {
  if(i %in% unique(sce$Subject)) {
    sce_m <- which(sce$Subject == i)
    meta_m <- which(subject_metadata_xlsx$Subject == i)
    sce$PCGA02_ParticipantID[sce_m] <- subject_metadata_xlsx$PCGA02_ParticipantID[meta_m]
    sce$Age_Baseline[sce_m] <- subject_metadata_xlsx$Age_Baseline_DSB[meta_m]
    sce$Cohort[sce_m] <- subject_metadata_xlsx$Cohort[meta_m]
    sce$PackYears[sce_m] <- subject_metadata_xlsx$PackYears[meta_m]
    sce$Ethnicity[sce_m] <- subject_metadata_xlsx$Ethnicity[meta_m]
    sce$Race[sce_m] <- subject_metadata_xlsx$Race[meta_m]
  }
}

## Sample Level Metadata
sample_metadata_xlsx <- metadata_xlsx[,c("PCGA02_ParentID", "Sample", "SampleType", "TimePoint", "Anatomic_Site", "Histology")]
sample_metadata_xlsx <- unique(sample_metadata_xlsx)
sample_metadata_xlsx$Anatomic_Site <- gsub(" ","", sample_metadata_xlsx$Anatomic_Site)
sce_full$Sample[which(sce_full$PCGA02_CellID_Template == "PCGA02_20152_3100054_1")] <- "54"
sce_full$TimePoint[which(sce_full$PCGA02_CellID_Template == "PCGA02_20152_3100054_1")] <- "T4" # Placeholder, as the other plate has this TP. Both plate's TP will be replaced by recent metadata
sce_full$PCGA02_ParentID <- ""
sce_full$SampleType <- "Biopsy"
sce_full$Anatomic_Site[which(sce_full$Anatomic_Site == "Trachea")] <- "TR"
coldata_sce_full <- data.frame(colData(sce_full))
coldata_sce_full_sample <- coldata_sce_full[,c("PCGA02_ParentID", "Sample", "SampleType", "TimePoint", "Anatomic_Site", "HistologyFull")]
coldata_sce_full_sample <- unique(coldata_sce_full_sample)
names(coldata_sce_full_sample)[6] <- "Histology"
coldata_sce_full_sample$TimePoint <- gsub("T","", coldata_sce_full_sample$TimePoint)
sce$Sample[which(sce$PCGA02_CellID_Template == "PCGA02_20152_3100054_1")] <- "54"
sce$TimePoint[which(sce$PCGA02_CellID_Template == "PCGA02_20152_3100054_1")] <- "T4" # Placeholder, as the other plate has this TP. Both plate's TP will be replaced by recent metadata
sce$PCGA02_ParentID <- ""
sce$SampleType <- "Biopsy"
sce$Anatomic_Site[which(sce$Anatomic_Site == "Trachea")] <- "TR"
coldata_sce <- data.frame(colData(sce))
coldata_sce_sample <- coldata_sce[,c("PCGA02_ParentID", "Sample", "SampleType", "TimePoint", "Anatomic_Site", "HistologyFull")]
coldata_sce_sample <- unique(coldata_sce_sample)
names(coldata_sce_sample)[6] <- "Histology"
coldata_sce_sample$TimePoint <- gsub("T","", coldata_sce_sample$TimePoint)

full_sample_reconciliation <- list()
for(i in sample_metadata_xlsx$Sample) {
  if(i %in% unique(coldata_sce_full_sample$Sample)) {
    xlsx_data <- sample_metadata_xlsx[which(sample_metadata_xlsx$Sample == i),]
    sce_full_data <- coldata_sce_full_sample[which(coldata_sce_full_sample$Sample == i),]
    full_sample_reconciliation[[i]] <- rbind(xlsx_data, sce_full_data)
  }
}
full_sample_reconciliation <- full_sample_reconciliation[!sapply(full_sample_reconciliation,is.null)]

# No discrepancies or missing values for SampleType. 
# Histology is labeled differently, even when the histologies are the same (e.g., due to capitalization, etc.). Changing those manaully.
table(sapply(full_sample_reconciliation, function(x) {length(unique(x[["PCGA02_ParentID"]]))}))
table(sapply(full_sample_reconciliation, function(x) {length(unique(x[["SampleType"]]))}))
table(sapply(full_sample_reconciliation, function(x) {length(unique(x[["TimePoint"]]))}))
table(sapply(full_sample_reconciliation, function(x) {length(unique(x[["Anatomic_Site"]]))}))
table(sapply(full_sample_reconciliation, function(x) {length(unique(x[["Histology"]]))}))
for(i in sample_metadata_xlsx$Sample) {
  if(i %in% unique(sce_full$Sample)) {
    sce_m <- which(sce_full$Sample == i)
    meta_m <- which(sample_metadata_xlsx$Sample == i)
    sce_full$PCGA02_ParentID[sce_m] <- sample_metadata_xlsx$PCGA02_ParentID[meta_m]
    sce_full$Anatomic_Site[sce_m] <- sample_metadata_xlsx$Anatomic_Site[meta_m]
    sce_full$TimePoint[sce_m] <- sample_metadata_xlsx$TimePoint[meta_m]
  }
}
sce_full$Anatomic_Site[which(sce_full$Anatomic_Site == "Unknown")] <- NA
sce_full$TimePoint <- paste0("T", sce_full$TimePoint)
sce_full$Histology[which(sce_full$Sample == "1630")] <- "Mild Dysplasia"
sce_full$HistologyFull[which(sce_full$Sample == "1630")] <- "Mild Dysplasia"
sce_full$Histology[which(sce_full$Sample == "1631")] <- "Mild Dysplasia"
sce_full$HistologyFull[which(sce_full$Sample == "1631")] <- "Mild Dysplasia"

nonimmune_sample_reconciliation <- list()
for(i in sample_metadata_xlsx$Sample) {
  if(i %in% unique(coldata_sce_sample$Sample)) {
    xlsx_data <- sample_metadata_xlsx[which(sample_metadata_xlsx$Sample == i),]
    sce_full_data <- coldata_sce_sample[which(coldata_sce_sample$Sample == i),]
    nonimmune_sample_reconciliation[[i]] <- rbind(xlsx_data, sce_full_data)
  }
}
nonimmune_sample_reconciliation <- nonimmune_sample_reconciliation[!sapply(nonimmune_sample_reconciliation,is.null)]

# No discrepancies or missing values for SampleType. 
# Histology is labeled differently, even when the histologies are the same (e.g., due to capitalization, etc.). Changing those manaully.
table(sapply(nonimmune_sample_reconciliation, function(x) {length(unique(x[["PCGA02_ParentID"]]))}))
table(sapply(nonimmune_sample_reconciliation, function(x) {length(unique(x[["SampleType"]]))}))
table(sapply(nonimmune_sample_reconciliation, function(x) {length(unique(x[["TimePoint"]]))}))
table(sapply(nonimmune_sample_reconciliation, function(x) {length(unique(x[["Anatomic_Site"]]))}))
table(sapply(nonimmune_sample_reconciliation, function(x) {length(unique(x[["Histology"]]))}))
for(i in sample_metadata_xlsx$Sample) {
  if(i %in% unique(sce$Sample)) {
    sce_m <- which(sce$Sample == i)
    meta_m <- which(sample_metadata_xlsx$Sample == i)
    sce$PCGA02_ParentID[sce_m] <- sample_metadata_xlsx$PCGA02_ParentID[meta_m]
    sce$Anatomic_Site[sce_m] <- sample_metadata_xlsx$Anatomic_Site[meta_m]
    sce$TimePoint[sce_m] <- sample_metadata_xlsx$TimePoint[meta_m]
  }
}
sce$Anatomic_Site[which(sce$Anatomic_Site == "Unknown")] <- NA
sce$TimePoint <- paste0("T", sce$TimePoint)
sce$Histology[which(sce$Sample == "1630")] <- "Mild Dysplasia"
sce$HistologyFull[which(sce$Sample == "1630")] <- "Mild Dysplasia"
sce$Histology[which(sce$Sample == "1631")] <- "Mild Dysplasia"
sce$HistologyFull[which(sce$Sample == "1631")] <- "Mild Dysplasia"

## Plate Level Metadata
plate_metadata_xlsx <- metadata_xlsx[,c("PCGA02_BiospecimenID", "experiment", "PCGA02_CellID_Template", "Batch", "R1", "R2", "Flow", "Dissociation_Enzyme", "Storage_Buffer")]
plate_metadata_xlsx <- unique(plate_metadata_xlsx)
plate_metadata_xlsx$Flow <- plate_metadata_xlsx$Dissociation_Enzyme
plate_metadata_xlsx$Dissociation_Enzyme <- NA
names(colData(sce_full))[26] <- "PCGA02_BiospecimenID"
sce_full$PCGA02_CellID_Template[which(sce_full$PCGA02_CellID_Template == "PCGA02-20152-3100079_2")] <- "PCGA02_20152_3100079_2"
sce_full$experiment[which(sce_full$experiment == "PCGA02-20152-3100079_1")] <- "PCGA02-20152-3000079_1"
sce_full$PCGA02_CellID_Template[which(sce_full$PCGA02_CellID_Template == "PCGA02_20152_3100079_1")] <- "PCGA02_20152_3000079_1"
sce_full$PCGA02_BiospecimenID[which(sce_full$PCGA02_CellID_Template == "PCGA02_20152_3000079_1")] <- "PCGA02_20152_3000079"
coldata_sce_full <- data.frame(colData(sce_full))
coldata_sce_full_plate <- coldata_sce_full[,c("PCGA02_BiospecimenID", "experiment", "PCGA02_CellID_Template", "Batch", "R1", "R2", "Flow", "Dissociation_Enzyme", "Storage_Buffer")]
coldata_sce_full_plate <- unique(coldata_sce_full_plate)
names(colData(sce))[26] <- "PCGA02_BiospecimenID"
sce$PCGA02_CellID_Template[which(sce$PCGA02_CellID_Template == "PCGA02-20152-3100079_2")] <- "PCGA02_20152_3100079_2"
sce$experiment[which(sce$experiment == "PCGA02-20152-3100079_1")] <- "PCGA02-20152-3000079_1"
sce$PCGA02_CellID_Template[which(sce$PCGA02_CellID_Template == "PCGA02_20152_3100079_1")] <- "PCGA02_20152_3000079_1"
sce$PCGA02_BiospecimenID[which(sce$PCGA02_CellID_Template == "PCGA02_20152_3000079_1")] <- "PCGA02_20152_3000079"
coldata_sce <- data.frame(colData(sce))
coldata_sce_plate <- coldata_sce[,c("PCGA02_BiospecimenID", "experiment", "PCGA02_CellID_Template", "Batch", "R1", "R2", "Flow", "Dissociation_Enzyme", "Storage_Buffer")]
coldata_sce_plate <- unique(coldata_sce_plate)

full_plate_reconciliation <- list()
for(i in plate_metadata_xlsx$PCGA02_CellID_Template) {
  if(i %in% unique(coldata_sce_full_plate$PCGA02_CellID_Template)) {
    xlsx_data <- plate_metadata_xlsx[which(plate_metadata_xlsx$PCGA02_CellID_Template == i),]
    sce_full_data <- coldata_sce_full_plate[which(coldata_sce_full_plate$PCGA02_CellID_Template == i),]
    full_plate_reconciliation[[i]] <- rbind(xlsx_data, sce_full_data)
  }
}
full_plate_reconciliation <- full_plate_reconciliation[!sapply(full_plate_reconciliation,is.null)]

# No discrepancies or missing values for PCGA02_BiospecimenID, experiment, Batch, R1, R2. 
# Dissociation_Enzyme not provided in Erin's spreadsheet, retained from prior versions.
# Dissociation_Enzyme more up to date in SCE object, so didn't overwrite it here.
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["PCGA02_BiospecimenID"]]))}))
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["experiment"]]))}))
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["PCGA02_CellID_Template"]]))}))
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["Batch"]]))}))
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["R1"]]))}))
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["R2"]]))}))
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["Flow"]]))}))
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["Dissociation_Enzyme"]]))}))
table(sapply(full_plate_reconciliation, function(x) {length(unique(x[["Storage_Buffer"]]))}))
for(i in plate_metadata_xlsx$PCGA02_CellID_Template) {
  if(i %in% unique(sce_full$PCGA02_CellID_Template)) {
    sce_m <- which(sce_full$PCGA02_CellID_Template == i)
    meta_m <- which(plate_metadata_xlsx$PCGA02_CellID_Template == i)
    sce_full$Flow[sce_m] <- plate_metadata_xlsx$Flow[meta_m]
  }
}

nonimmune_plate_reconciliation <- list()
for(i in plate_metadata_xlsx$PCGA02_CellID_Template) {
  if(i %in% unique(coldata_sce_plate$PCGA02_CellID_Template)) {
    xlsx_data <- plate_metadata_xlsx[which(plate_metadata_xlsx$PCGA02_CellID_Template == i),]
    sce_full_data <- coldata_sce_plate[which(coldata_sce_plate$PCGA02_CellID_Template == i),]
    nonimmune_plate_reconciliation[[i]] <- rbind(xlsx_data, sce_full_data)
  }
}
nonimmune_plate_reconciliation <- nonimmune_plate_reconciliation[!sapply(nonimmune_plate_reconciliation,is.null)]

# No discrepancies or missing values for PCGA02_BiospecimenID, experiment, Batch, R1, R2. 
# Dissociation_Enzyme not provided in Erin's spreadsheet, retained from prior versions.
# Dissociation_Enzyme more up to date in SCE object, so didn't overwrite it here.
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["PCGA02_BiospecimenID"]]))}))
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["experiment"]]))}))
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["PCGA02_CellID_Template"]]))}))
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["Batch"]]))}))
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["R1"]]))}))
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["R2"]]))}))
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["Flow"]]))}))
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["Dissociation_Enzyme"]]))}))
table(sapply(nonimmune_plate_reconciliation, function(x) {length(unique(x[["Storage_Buffer"]]))}))
for(i in plate_metadata_xlsx$PCGA02_CellID_Template) {
  if(i %in% unique(sce$PCGA02_CellID_Template)) {
    sce_m <- which(sce$PCGA02_CellID_Template == i)
    meta_m <- which(plate_metadata_xlsx$PCGA02_CellID_Template == i)
    sce$Flow[sce_m] <- plate_metadata_xlsx$Flow[meta_m]
  }
}

altExp(sce_full)$PCGA02_ParticipantID <- sce_full$PCGA02_ParticipantID
sce_full$Age_Baseline <- as.numeric(sce_full$Age_Baseline)
altExp(sce_full)$Age_Baseline <- sce_full$Age_Baseline
altExp(sce_full)$PackYears <- sce_full$PackYears
altExp(sce_full)$Ethnicity <- sce_full$Ethnicity
altExp(sce_full)$Race <- sce_full$Race
altExp(sce_full)$PCGA02_ParentID <- sce_full$PCGA02_ParentID
altExp(sce_full)$SampleType <- sce_full$SampleType
altExp(sce_full)$TimePoint <- sce_full$TimePoint
altExp(sce_full)$Anatomic_Site <- sce_full$Anatomic_Site
altExp(sce_full)$Histology <- sce_full$Histology
altExp(sce_full)$HistologyFull <- sce_full$HistologyFull
altExp(sce_full)$PCGA02_BiospecimenID <- sce_full$PCGA02_BiospecimenID
altExp(sce_full)$experiment <- sce_full$experiment
altExp(sce_full)$PCGA02_CellID_Template <- sce_full$PCGA02_CellID_Template
altExp(sce_full)$Flow <- sce_full$Flow
names(colData(sce_full))[44] <- "PCGA02_BiospecimenID_Batch"
sce_full$PCGA02_BiospecimenID_Batch <- paste(sce_full$PCGA02_BiospecimenID, "batch", sce_full$Batch, sep = "_")
altExp(sce_full)$PCGA02_BiospecimenID_Batch <- sce_full$PCGA02_BiospecimenID_Batch
sce_full$Subject_TimePoint <- paste(sce_full$Subject, sce_full$TimePoint, sep = "_")
altExp(sce_full)$Subject_TimePoint <- sce_full$Subject_TimePoint
altExp(sce_full)$PCGA02_LongID <- c()
sce_full$Batch <- as.character(as.numeric(sce_full$Batch))
altExp(sce_full)$Batch <- sce_full$Batch
altExp(sce_full)$PCGA02_LongID_Batch <- c()
altExp(sce_full)$Smoking_Status <- sce_full$Smoking_Status
altExp(sce_full)$Sex <- sce_full$Sex
altExp(sce_full)$celda_sample_label <- altExp(sce_full)$PCGA02_BiospecimenID
altExp(sce_full)$colnames <- colnames(sce_full)
altExp(sce_full)$HGBasalPercent <- sce_full$HGBasalPercent
altExp(sce_full)$HGBasalPercentMixed <- sce_full$HGBasalPercentMixed
altExp(sce_full)$TotalBasalNumbers <- sce_full$TotalBasalNumbers
altExp(sce_full)$TotalBasalNumbersMixed <- sce_full$TotalBasalNumbersMixed
altExp(sce_full)$BasalGroup <- sce_full$BasalGroup
altExp(sce_full)$BasalGroupMixed <- sce_full$BasalGroupMixed
altExp(sce_full)$CellType_Cluster <- sce_full$CellType_Cluster

altExp(sce)$PCGA02_ParticipantID <- sce$PCGA02_ParticipantID
sce$Age_Baseline <- as.numeric(sce$Age_Baseline)
altExp(sce)$Age_Baseline <- sce$Age_Baseline
altExp(sce)$PackYears <- sce$PackYears
altExp(sce)$Ethnicity <- sce$Ethnicity
altExp(sce)$Race <- sce$Race
altExp(sce)$PCGA02_ParentID <- sce$PCGA02_ParentID
altExp(sce)$SampleType <- sce$SampleType
altExp(sce)$TimePoint <- sce$TimePoint
altExp(sce)$Anatomic_Site <- sce$Anatomic_Site
altExp(sce)$Histology <- sce$Histology
altExp(sce)$HistologyFull <- sce$HistologyFull
altExp(sce)$PCGA02_BiospecimenID <- sce$PCGA02_BiospecimenID
altExp(sce)$experiment <- sce$experiment
altExp(sce)$PCGA02_CellID_Template <- sce$PCGA02_CellID_Template
altExp(sce)$Flow <- sce$Flow
sce$PCGA02_BiospecimenID_Batch <- paste(sce$PCGA02_BiospecimenID, "batch", sce$Batch, sep = "_")
altExp(sce)$PCGA02_BiospecimenID_Batch <- sce$PCGA02_BiospecimenID_Batch
sce$Subject_TimePoint <- paste(sce$Subject, sce$TimePoint, sep = "_")
altExp(sce)$Subject_TimePoint <- sce$Subject_TimePoint
altExp(sce)$PCGA02_LongID <- c()
sce$Batch <- as.character(as.numeric(sce$Batch))
altExp(sce)$Batch <- as.character(as.numeric(altExp(sce)$Batch))
altExp(sce)$PCGA02_LongID_Batch <- c()
altExp(sce)$Smoking_Status <- sce$Smoking_Status
altExp(sce)$Sex <- sce$Sex
altExp(sce)$celda_sample_label <- altExp(sce)$PCGA02_BiospecimenID
altExp(sce)$colnames <- colnames(sce)
altExp(sce)$KCluster <- sce$KCluster
altExp(sce)$Cluster_CellType <- sce$Cluster_CellType
altExp(sce)$celda_cell_cluster_allcells <- sce$celda_cell_cluster_allcells
altExp(sce)$CellType_Cluster <- sce$CellType_Cluster
sce$HistologyCondense <- c()
altExp(sce)$BasalCluster <- sce$BasalCluster
altExp(sce)$BasalClusterMixed <- sce$BasalClusterMixed
altExp(sce)$HGBasalPercent <- sce$HGBasalPercent
altExp(sce)$HGBasalPercentMixed <- sce$HGBasalPercentMixed
altExp(sce)$TotalBasalNumbers <- sce$TotalBasalNumbers
altExp(sce)$TotalBasalNumbersMixed <- sce$TotalBasalNumbersMixed
altExp(sce)$BasalGroup <- sce$BasalGroup
altExp(sce)$BasalGroupMixed <- sce$BasalGroupMixed

### 02/19/25: This will only happen after some additional metadata reconciliation
### 02/29/25: The objects below are currently saved with the suffix "_PostMetadataRec_020325"
saveRDS(sce_full, "121224_AllCells_SCE_BasalGroup.rds")
saveRDS(sce, "121224_Nonimmune_SCE_BasalGroup.rds")

save.image(paste0(format(Sys.Date(),"%m%d%y"),"_Metadata_Rec.RData"))
