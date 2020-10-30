#This script details the pre-processing step done on scRNA-seq before guided clustering was performed.
#The scRNA-seq data used is publicly available from Tucker et. al. which can be downloaded [here](https://singlecell.broadinstitute.org/single_cell/study/SCP498/transcriptional-and-cellular-diversity-of-the-human-heart#study-download).
library(Seurat)
library(Matrix)

#Read normalised 10X scRNA-seq dataset. The dataset contains 287,269 cells from the 4 chambers of normal human heart samples (LV, LA, RV, RA) and 33,694 genes.
sc.counts <- Read10X("/Volumes/Seagate Backup Plus Drive/UROPS/Datasets/Single-Cell Data")

#Sample selection step: For the purpose of the project, retain only cells from the LV.
LV.cells <- sc.counts[,substr(colnames(sc.counts), 1, 2) == "LV"]

#Pre-process genes to ensure we only cluster genes with meaningful data. Retain genes that meet two criteria: 
#1) Has counts in at least 10 cells and 2) Belongs to the top 30% highly expressed genes.
LV.cells <- LV.cells[rowSums(LV.cells != 0) >= 10, ,drop = FALSE]

#Due to RAM limitations, only retain the top 30% highly-expressed genes.  
LV.cells <- LV.cells[rev(order(rowSums(LV.cells))),]
LV.cells <- LV.cells[1:7494,]
