#This script details the pre-processing steps done on bulk RNA-Seq and unnormalised scRNA-seq before bulk tissue deconvolution was performed using the Bisque R package. 
#Pre-processing was done due to RAM limitations.

library(Seurat)
library(Biobase)
library(data.table)

#Load unnormalised scRNA-seq data and extract counts.
sc.data <-  ReadH5AD("path/healthy_human_4chamber_map_unnormalized_V3.h5ad")
sc.counts <- sc.data@assays$RNA@counts

#Load gene name to gene ID matrix. Replace gene names with gene IDs in scRNA-seq to make it comparable with gene IDs in bulk RNA-seq.
sc.genes <- read.table("path/genes.tsv")
rownames(sc.counts) <- sc.genes$V1

#Label cells in scRNA-seq by individual IDs for Bisque to perform deconvolution.
individual.ids <- sc.data@meta.data$biological.individual
individual.ids <- t(individual.ids)
sc.counts <- rbind(sc.counts, individual.ids)

#Load guided clustering output. Cells are annotated with cell-types for Bisque to reference. 
#Note that only a subset of the original output was used due to RAM limitations. The subset used contains 10% of the cells from each cell type, randomly selected via proportional sampling.
sampled.cells <- readRDS("path/SampledCells.RData")

#Subset scRNA-seq data to retain only sampled cells.
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
