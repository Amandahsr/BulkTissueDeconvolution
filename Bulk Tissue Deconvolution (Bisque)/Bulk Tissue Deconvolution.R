#This script details the bulk tissue deconvolution step using the Bisque Reference-Based Decomposition Model.
#Bisque requires 4 inputs for the Reference-Based Decompostion Model:
#1. Unnormalised bulk RNA-seq data in ExpressionSet format.
#2. Unnormalised scRNA-seq data in ExpressionSet format.
#3. Cell-type proportions from guided clustering step.
#4. DEGs from guided clustering step (OPTIONAL).

library(data.table)
library(Biobase)
library(BisqueRNA)

#Read unnormalised bulk RNA-seq dataset.
#Dataset contains 62 bulk tissue samples with 22,474 genes. There are 33 normal heart samples and 29 DCM heart samples.
bulk.data <- read.delim("/Volumes/Seagate Backup Plus Drive/Datasets/Bulk-Tissue/raw_count.txt")

#Gene IDs in bulk RNA-seq data are processed to remove decimal places for comparison with genes in scRNA-seq data.
rownames(bulk.data) <- gsub("\\.[0-9]","", rownames(bulk.data))

#Read bulk tissue sample info.
experiment.labels <- read.table("/Volumes/Seagate Backup Plus Drive/Datasets/Bulk-Tissue/grouping.txt")

#Retain only normal heart samples in bulk RNA-seq data.
control.labels <- experiment.labels[experiment.labels$V2 == 'NF',]
bulk.control <- bulk.data[,intersect(control.labels$V1, colnames(bulk.data))]

#Convert bulk RNA-seq data to ExpressionSet format for bulk tissue deconvolution.
bulk.control <- as.matrix(bulk.control)
bulk.exp <- ExpressionSet(bulk.control)

#Read DEGs from guided clustering to use as markers for bulk tissue deconvolution.
deg.genes <- read.table("/Volumes/Seagate Backup Plus Drive/Datasets/DEG Analysis/Unique_DEGs2.csv")
gene.markers <- subset(genes, V2 %in% deg.genes$V1)
gene.ids <- as.vector(gene.markers$V1)

#Perform bulk tissue deconvolution using the Bisque Reference-Based Decomposition Model.
#use.overlap was set to FALSE because there are no common samples between scRNA-seq and bulk RNA-seq data.
control.est <- BisqueRNA::ReferenceBasedDecomposition(bulk.exp, sc.exp, markers = gene.ids, use.overlap = FALSE)



