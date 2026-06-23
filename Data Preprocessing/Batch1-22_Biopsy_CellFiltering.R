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

setwd("/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/Analysis/PCGA_Biopsy_Paper/Batch1-22_Biopsies/")

sce_unf <- readRDS("032624_Biopsies_Batches1-22_Unfiltered.rds")
dim(sce_unf) # 60716 14784

# ERCC Correlation
control = read.table("/restricted/projectnb/pcga/Conor/references/ercc_conc.txt",
                     sep = '\t', header = T)
rownames(control) <- control$ERCC.ID
sce_unf$ERCCcorrelation <- apply(counts(sce_unf)[rownames(control),], 2, function(e){
  return(summary(lm(e ~ control$concentration.in.Mix.1..attomoles.ul.))$adj.r.squared)})
table(is.na(sce_unf$ERCCcorrelation)) # FALSE: 11320  TRUE: 3464
sce_unf$ERCCcorrelation[is.na(sce_unf$ERCCcorrelation)] <- 0
sce_unf$ERCCcounts <- apply(counts(sce_unf)[rownames(control),], 2, sum)

# Subset to protein coding genes
protein_coding <- which(rowData(sce_unf)$Gene_biotype == "protein_coding")
sce_unf <- sce_unf[protein_coding,]
dim(sce_unf) # 19941 14784
mt_nopseudo <- setdiff(seq_len(nrow(sce_unf)), grep("MTRNR", rowData(sce_unf)$Gene_name))
sce_unf <- sce_unf[mt_nopseudo,]
dim(sce_unf) # 19930 14784

# Mitochondrial Count Percent
mt_genes <- grep("^MT-",rowData(sce_unf)$Gene_name)
protein_genes <- grep("protein_coding", rowData(sce_unf)$Gene_biotype)
mt_genes <- intersect(mt_genes, protein_genes)
rowData(sce_unf)$Gene_name[mt_genes] # "MT-ND6"  "MT-CO2"  "MT-CYB"  "MT-ND2"  "MT-ND5"  "MT-CO1"  "MT-ND3"  "MT-ND4"  "MT-ND1"  "MT-ATP6" "MT-CO3"  "MT-ND4L" "MT-ATP8"
geneSet <- list("mt.percent_runPerCellQC" = rownames(sce_unf)[mt_genes])

# Filter empty wells
cd <- data.frame(colData(sce_unf))
empty_well <- cd[which(cd$number_of_cells %in% c("101","112","801","812")),]
empty_well_df <- empty_well %>% group_by(PCGA02_CellID_Template) %>% summarize(minwell = min(reads))
length(unique(empty_well$reads)) # 616
length(empty_well$reads) # 616
empty_wells <- which(cd$number_of_cells %in% c("101","112","801","812") & cd$reads %in% empty_well_df$minwell)
nonempty_wells <- setdiff(seq_len(ncol(sce_unf)), empty_wells)
sce_unf <- sce_unf[,nonempty_wells]
dim(sce_unf) # 19930 14630

# Filter genes with no counts
hascts <- which(rowSums(counts(sce_unf)) > 0)
sce_unf <- sce_unf[hascts,]
dim(sce_unf) # 18485 14630

# Run QC
sce_unf <- runPerCellQC(sce_unf, geneSetList = geneSet)
sce_unf <- runCellQC(inSCE = sce_unf, algorithms = c("scDblFinder","cxds","bcds","cxds_bcds_hybrid","doubletFinder"))
names(colData(sce_unf))[59] <- "mt_percent"
names(colData(sce_unf)) <- gsub("percent.top", "percent_top", names(colData(sce_unf)))
saveRDS(sce_unf, paste0(format(Sys.Date(),"%m%d%y"), "_Biopsies_Batches1-22_Unfiltered_QC+Doublets.rds"))

# Filter poor quality cells
filter_thresholds <- which(sce_unf$sum < 500 | sce_unf$mt_percent > 40 | sce_unf$percent_top_50 > 75)
sce <- sce_unf[,-filter_thresholds]
dim(sce) # 18485  9435
View(table(sce$PCGA02_CellID_Template))
# Remove plates with fewer than 10 cells
bad_plates <- c("PCGA02_10417_3002834_1",
                "PCGA02_20005_3100004_1",
                "PCGA02_20127_3100014_1",
                "PCGA02_20149_3100033_1",
                "PCGA02_10420_3001724_1",
                "PCGA02_20137_3137105_1",
                "PCGA02_10405_3001778_1")
sce <- sce[,which(!sce$PCGA02_CellID_Template %in% bad_plates)]
dim(sce) # 18485  9405

saveRDS(sce, paste0(format(Sys.Date(),"%m%d%y"), "_Biopsies_Batches1-22_PCGenes_NoMTPseudo_Filtered_UMI>500_MT<0.4_p50<0.75.rds"))

save.image("Batch1-22_CellFiltering.RData")

