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
library(mgsub)

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies")

sce_batch1.21 <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-21_Biopsies/030424_Biopsies_Batches1-21_Unfiltered_noChangedFlow.rds")

sce_batch22 <- readRDS("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/PCGA_PML-Biopsy_Batch22/Alignment/scruff_object_final_03082024.rds")
sce_batch22 <- sce_batch22[,-grep("HTA3", sce_batch22$experiment)] # remove HTAN Colorado plates
sce_batch22 <- sce_batch22[,-grep("3601839", sce_batch22$experiment)] # remove nasal brush plate
sce_batch22$PCGA02_LongID <- sce_batch22$experiment
sce_batch22$PCGA02_CellID_Template <- paste0(sce_batch22$PCGA02_LongID,"_1")
sce_batch22$Batch <- "22"
sce_batch22$Subject <- rep(c("420","405","203","405","163","163","417","417"), each = 96)
sce_batch22$Sample <- rep(c("3420","1778","1868","1778","63105","63106","3724","3724"), each = 96)
sce_batch22$Flow <- rep(c("CD45-","CD45-","CD45-","CD45+","CD45+","CD45-","CD45-","CD45+"), each = 96)
batch22_metadata <- read_excel("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Batch_Metadata/Batch22_Plate_Information.xlsx", col_names = TRUE)
m <- match(sce_batch22$experiment, batch22_metadata$experiment)
sce_batch22$R1 <- batch22_metadata$R1[m]
sce_batch22$R2 <- batch22_metadata$R2[m]
batch22_core <- read_excel("PCGA2.0_Batch22_Project File_11022023.xlsx", col_names = TRUE)
m <- match(sce_batch22$experiment, batch22_core$`Biospecimen ID`)
sce_batch22$PCGA02_site <- batch22_core$Site_Text[m]
sce_batch22$TimePoint <- batch22_core$Timepoint[m]
sce_batch22$Dissociation_Enzyme <- "Trypsin"
sce_batch22$Storage_Buffer <- batch22_core$Buffer[m]
sce_batch22$Age_Baseline <- rep(c("57","61","69","61","70","70","58","58"), each = 96)
sce_batch22$Sex <- rep(c("Female","Male","Male","Male","Male","Male","Male","Male"), each = 96)
sce_batch22$Smoking_Status <- rep(c("Current","Former","Former","Former","Former","Former","Current","Current"), each = 96)
sce_batch22$Anatomic_Site <- batch22_core$AnatomicSite_Text[m]
sce_batch22$Histology <- batch22_core$Histology[m]
sce_batch22$Cohort <- batch22_core$PCGA02_Cohort[m]
sce_batch22$PCGA02_LongID_Batch <- paste0(sce_batch22$PCGA02_LongID,"_22")
sce_batch22$Subject_TimePoint <- paste(sce_batch22$Subject, sce_batch22$TimePoint, sep = "_")
sce_batch22$PackYears <- rep(c(37, 43, 180, 43, 50, 50, 54, 54), each = 96)
sce_batch22$Ethnicity <- NA
sce_batch22$Race <- NA
colnames(sce_batch22) <- paste0(sce_batch22$PCGA02_CellID_Template, rep(sprintf("%02d",1:96), times = 8))

sce_batch1.22 <- cbind(sce_batch1.21, sce_batch22)
subject_metadata <- read_excel("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Batch_Metadata/scRNA_allBatches_compiled_metadata_2023-12-11_ConorEdit.xlsx", 
                               sheet = "Patient Demographics", col_names = TRUE)
m <- match(sce_batch1.22$Subject, subject_metadata$Local_Subject_ID)
for(i in seq_along(m)) {
  if(!is.na(m[i])) {
    sce_batch1.22$Age_Baseline[i] <- subject_metadata$Age_Baseline[m[i]]
    sce_batch1.22$Sex[i] <- subject_metadata$Sex[m[i]]
    sce_batch1.22$Smoking_Status[i] <- subject_metadata$Smoker[m[i]]
    sce_batch1.22$PackYears[i] <- subject_metadata$PackYears[m[i]]
  } else {
    sce_batch1.22$Age_Baseline[i] <- sce_batch1.22$Age_Baseline[i]
    sce_batch1.22$Sex[i] <- sce_batch1.22$Sex[i]
    sce_batch1.22$Smoking_Status[i] <- sce_batch1.22$Smoking_Status[i]
    sce_batch1.22$PackYears[i] <- sce_batch1.22$PackYears[i]
  }
}

sce_batch1.22$Batch <- factor(sce_batch1.22$Batch, levels = c("1" ,"2","3","4","5","7","8","9","10","11","12","19","20","21","22"), ordered = TRUE)
sce_batch1.22$Sex <- factor(sce_batch1.22$Sex, levels = c("Female","Male"), ordered = TRUE)
sce_batch1.22$Smoking_Status <- factor(sce_batch1.22$Smoking_Status, levels = c("Former","Current"), ordered = TRUE)
sce_batch1.22$Histology <- factor(sce_batch1.22$Histology, levels = c("Normal","Denuded Bronchial Mucosa","Inflammation","Basal Cell Hyperplasia","Metaplasia", 
                                                        "Mild Dysplasia","Moderate Dysplasia","Severe Dysplasia","CIS","LUSC"), ordered = TRUE)
sce_batch1.22$number_of_cells <- as.character(sce_batch1.22$number_of_cells)
names(rowData(sce_batch1.22)) <- c("Gene_id","Gene_name","Gene_biotype","Seqid","Start","End","Strand","Source")

saveRDS(sce_batch1.22, paste0(format(Sys.Date(),"%m%d%y"), "_Biopsies_Batches1-22_Unfiltered.rds"))

save.image("Batch1-22_DataManagement.RData")

