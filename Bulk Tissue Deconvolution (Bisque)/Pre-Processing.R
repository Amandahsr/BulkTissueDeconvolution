#This script details the pre-processing step done on scRNA-seq before bulk tissue deconvolution was performed.

library(data.table)
library(Seurat)
library(Biobase)

#Read unnormalised scRNA-seq dataset.
#Dataset contains 287,269 cells from the 4 chambers of normal human heart samples (LV, LA, RV, RA) and 33,694 genes.
sc.data <-  ReadH5AD("/Volumes/Seagate Backup Plus Drive/UROPS/Datasets/Single-Cell Data/healthy_human_4chamber_map_unnormalized_V3.h5ad")
sc.counts <- sc.data@assays$RNA@counts

#Replace gene names with gene IDs in scRNA-seq because genes are labelled by IDs in bulk RNA-seq.
sc.genes <- read.table("/Volumes/Seagate Backup Plus Drive/UROPS/Datasets/Single-Cell Data/genes.tsv")
rownames(sc.counts) <- sc.genes$V1

#Label cells by individual IDs.
individual.ids <- sc.data@meta.data$biological.individual
individual.ids <- t(individual.ids)
sc.counts <- rbind(sc.counts, individual.ids)

#Read data for LV cells that were randomly sampled using proportional sampling.
#Contains 10% of the cells from each cell type.
#Cells are labelled by cell type.
sampled.cells <- read.csv("/Volumes/Seagate Backup Plus Drive/UROPS/Datasets/Bulk Tissue Data/Proportional_sampling_cells.csv")

#Subset scRNA-seq data to only sampled cells.
sc.counts <- sc.counts[,match(sampled.cells$SampleID,colnames(sc.counts))]

#Convert scRNA-seq data to ExpressionSet format for bulk tissue deconvolution.
sc.matrix.counts <- as.matrix(sc.counts[1:33694,])
individual.ids <- data.frame(sc.counts[33695,])
rownames(individual.ids) <- NULL
colnames(individual.ids) <- 'SubjectName'
sample.ids <- colnames(sc.counts)
sc.pheno <- data.frame(check.names=F, check.rows=F, stringsAsFactors=F,row.names=sample.ids,SubjectName = individual.ids, cellType = sampled.cells$cellType)
sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"), row.names=c("SubjectName", "cellType"))
sc.pdata <- new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta)
sc.exp <- ExpressionSet(assayData = sc.matrix.counts, phenoData = sc.pdata)
