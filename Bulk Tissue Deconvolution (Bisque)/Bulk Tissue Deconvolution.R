#This script details the pre-processing of Bulk RNbulk tissue deconvolution step using the Bisque Reference-Based Decomposition Model.
#Bisque requires 4 inputs: 1) Unnormalised bulk RNA-seq data, 2) Unnormalised scRNA-seq data, 3) Cell-type proportions from guided clustering results, 4) DEGs from guided clustering step (OPTIONAL).library(BisqueRNA)
library(Bisque)
library(Biobase)
library(data.table)

#Load the required 4 inputs. The ExpressionSet of scRNA-seq would already contain the cell-type proportions from guided clustering results.
sc.exp <- readRDS("path/scEXP.RData")
bulk.exp <- readRDS("path/bulkEXP.RData")
deg.genes <- readRDS("path/ClusteringDEG.RData")

#Load gene name to gene ID matrix. Replace DEG names with gene IDs to make them comparable with gene IDs in scRNA-seq and bulk RNA-seq data.
sc.genes <- readRDS("path/GeneNametoID.RData")
gene.markers <- subset(sc.genes, V2 %in% deg.genes$V1)
gene.ids <- as.vector(gene.markers$V1)

#Perform bulk tissue deconvolution using the Bisque Reference-Based Decomposition Model. use.overlap was set to FALSE because there are no common samples between scRNA-seq and bulk RNA-seq data.
control.est <- BisqueRNA::ReferenceBasedDecomposition(bulk.exp, sc.exp, markers = gene.ids, use.overlap = FALSE)



