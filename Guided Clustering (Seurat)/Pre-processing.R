#This script details the pre-processing step done on scRNA-seq before guided clustering was performed.

library(Seurat)
library(Matrix)

#Read normalised scRNA-seq dataset.
#Dataset contains 287,269 cells from the 4 chambers of normal human heart samples (LV, LA, RV, RA) and 33,694 genes.
sc.counts <- Read10X("/Volumes/Seagate Backup Plus Drive/UROPS/Datasets/Single-Cell Data")

#Sample selection step:
#Retain only cells from the LV for guided clustering.
LV.cells <- sc.counts[,substr(colnames(sc.counts), 1, 2) == "LV"]

#Genes selection step:
#Retain genes that meet two criteria: 
#1) Has counts in at least 10 cells and 2) Belongs to the top 30% highly expressed genes
LV.cells <- LV.cells[rowSums(LV.cells != 0) >= 10, ,drop = FALSE]
LV.cells <- LV.cells[rev(order(rowSums(LV.cells))),]
LV.cells <- LV.cells[1:7494,]
