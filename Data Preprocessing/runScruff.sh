#!/bin/bash -l

#$ -S /bin/bash

#$ -P pcga

#$ -j y

#$ -l h_rt=36:00:00

#$ -l mem_per_core=8G

#$ -pe omp 8

echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="

# Find R script dynamically using find
start=`date +%s`
code=$(find -O3 /restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/ -name "reportScruffQC.R")
end=`date +%s`
DIFF="The runtime for the 'find' command was "$(( $end - $start ))" seconds."
echo $DIFF

out_dir="/restricted/projectnb/pcga/PCGA2_BronchialPML_SingleCell_Biopsies+Brushes/PCGA_PML-Biopsy_Batch${BATCH}/Alignment/"
mkdir -p $out_dir
cd $out_dir
source ~/.bashrc
module load R/4.0.2
module load pandoc/2.5

Rscript $code \
-P $out_dir \
-B ${BATCH} \
-F "${FILE}"
