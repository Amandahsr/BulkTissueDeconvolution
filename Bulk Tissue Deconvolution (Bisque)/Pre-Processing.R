#This script details the pre-processing steps done on bulk RNA-seq and unnormalised scRNA-seq before bulk tissue deconvolution was performed using the Bisque R package. 
#Pre-processing is required because guided clustering was performed using a subset of the original scRNA-seq data, and the project is focused on deconvoluting normal heart tissue samples
library(Seurat)
library(Biobase)
library(data.table)

#Pre-processing of scRNA-seq data
############################################################################
#Load unnormalised scRNA-seq data and extract counts.
sc.data <-  ReadH5AD("path/healthy_human_4chamber_map_unnormalized_V3.h5ad")
sc.counts <- sc.data@assays$RNA@counts

#Load gene name to gene ID matrix. Replace gene names with gene IDs in scRNA-seq to make it comparable with gene IDs in bulk RNA-seq.
sc.genes <- readRDS("path/GeneNametoID.RData")
rownames(sc.counts) <- sc.genes$V1

#Label cells in scRNA-seq by individual IDs for Bisque to perform deconvolution.
individual.ids <- sc.data@meta.data$biological.individual
individual.ids <- t(individual.ids)
sc.counts <- rbind(sc.counts, individual.ids)

#Load guided clustering output. Cells are annotated with cell-types for Bisque to reference. 
#Note that only a subset of the original output was used due to RAM limitations. 
#The subset used contains 10% of the cells from each cell type, randomly selected via proportional sampling.
sampled.cells <- readRDS("path/SampledCells.RData")

#Subset scRNA-seq data to retain only sampled cells.
sc.counts <- sc.counts[,match(sampled.cells$SampleID,colnames(sc.counts))]

#Convert scRNA-seq data to ExpressionSet format for bulk tissue deconvolution. Incorporate guided clustering results into the ExpressionSet by annotating cells with cell-types.
sc.matrix.counts <- as.matrix(sc.counts[1:33694,])
individual.ids <- data.frame(sc.counts[33695,])
rownames(individual.ids) <- NULL
colnames(individual.ids) <- 'SubjectName'
sample.ids <- colnames(sc.counts)
sc.pheno <- data.frame(check.names=F, check.rows=F, stringsAsFactors=F,row.names=sample.ids,SubjectName = individual.ids, cellType = sampled.cells$cellType)
sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"), row.names=c("SubjectName", "cellType"))
sc.pdata <- new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta)
sc.exp <- ExpressionSet(assayData = sc.matrix.counts, phenoData = sc.pdata)

#Pre-processing of bulk RNA-seq data
############################################################################
#Load bulk RNA-seq data.
bulk.data <- read.delim("path/raw_count.txt")

#Decimal places of gene IDs in bulk RNA-seq data are removed for comparison with gene IDs in scRNA-seq data.
rownames(bulk.data) <- gsub("\\.[0-9]","", rownames(bulk.data))

#For the purposes of the project, read bulk tissue sample info and retain only normal heart samples in bulk RNA-seq data.
experiment.labels <- read.table("/Volumes/Seagate Backup Plus Drive/UROPS/Datasets/Bulk Tissue Data/grouping.txt")
control.labels <- experiment.labels[experiment.labels$V2 == 'NF',]
bulk.control <- bulk.data[,intersect(control.labels$V1, colnames(bulk.data))]

#Convert bulk RNA-seq data to ExpressionSet format for bulk tissue deconvolution.
bulk.control <- as.matrix(bulk.control)
bulk.exp <- ExpressionSet(bulk.control)
