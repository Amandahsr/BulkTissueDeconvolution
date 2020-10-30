#This script details the full steps of bulk tissue deconvolution on bulk RNA-seq data.
library(Seurat)
library(BisqueRNA)
library(Biobase)
library(data.table)

##1.Pre-processing steps:
##############################################################################################################
#Pre-processing of scRNA-seq data
#################################
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
####################################
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

##2.Bulk Tissue Deconvolution steps:
##############################################################################################################
#Load gene name to gene ID matrix. Replace DEG names with gene IDs to make them comparable with gene IDs in scRNA-seq and bulk RNA-seq data.
sc.genes <- readRDS("path/GeneNametoID.RData")
gene.markers <- subset(sc.genes, V2 %in% deg.genes$V1)
gene.ids <- as.vector(gene.markers$V1)

#Perform bulk tissue deconvolution using the Bisque Reference-Based Decomposition Model. use.overlap was set to FALSE because there are no common samples between scRNA-seq and bulk RNA-seq data.
control.est <- BisqueRNA::ReferenceBasedDecomposition(bulk.exp, sc.exp, markers = gene.ids, use.overlap = FALSE)

##3.Output steps:
##############################################################################################################
#Extract cell-type proportions from each sample.
control.prop <- data.frame(control.est$bulk.props)

#Calculate mean of cell-type proportions across samples. 
#Mean cell-type proportions are used as comparison with stained normal heart tissue to evaluate the accuracy of the Bisque algorithm.
control.avgprop <- data.frame(rowMeans(control.prop))
