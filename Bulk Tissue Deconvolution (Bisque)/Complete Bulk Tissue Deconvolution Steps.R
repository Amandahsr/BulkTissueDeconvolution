#This script details the full steps of bulk tissue deconvolution on bulk RNA-seq data.

library(data.table)
library(Seurat)
library(Biobase)
library(BisqueRNA)

##1.Pre-processing steps:
#######################################################
#Read unnormalised scRNA-seq dataset.
#Dataset contains 287,269 cells from the 4 chambers of normal human heart samples (LV, LA, RV, RA) and 33,694 genes.
sc.data <-  ReadH5AD("/Volumes/Seagate Backup Plus Drive/Datasets/Single-Cell Data/healthy_human_4chamber_map_unnormalized_V3.h5ad")
sc.counts <- sc.data@assays$RNA@counts

#Replace gene names with gene IDs in scRNA-seq because genes are labelled by IDs in bulk RNA-seq.
sc.genes <- read.table("/Volumes/Seagate Backup Plus Drive/Datasets/Single-Cell Data/genes.tsv")
rownames(sc.counts) <- sc.genes$V1

#Label cells by individual IDs.
individual.ids <- sc.data@meta.data$biological.individual
individual.ids <- t(individual.ids)
sc.counts <- rbind(sc.counts, individual.ids)

#Read data for LV cells that were randomly sampled using proportional sampling.
#Contains 10% of the cells from each cell type.
#Cells are labelled by cell type.
sampled.cells <- read.csv("/Volumes/Seagate Backup Plus Drive/Datasets/Bulk-Tissue/Proportional_sampling_cells.csv")

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

##2.Bulk Tissue Deconvolution steps:
#######################################################
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

##3.Output steps:
#######################################################
#Extract cell-type proportions in each sample.
control.prop <- data.frame(control.est$bulk.props)

#Calculate mean of cell-type proportions across samples.
#Mean cell-type proportions are used as comparison with stained normal heart tissue to evaluate the accuracy of the Bisque algorithm.
control.avgprop <- data.frame(rowMeans(control.prop))
