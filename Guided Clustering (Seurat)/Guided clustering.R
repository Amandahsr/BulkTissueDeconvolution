#This script details the guided clustering step on pre-processed scRNA-seq data.

library(Seurat)

#Initialise UMI counts as Seurat object for clustering. 
#min.cells and min.features are set to 0 because pre-processing was already performed.
seurat.data <- CreateSeuratObject(LV.cells, project = 'Seurat Guided Clustering', min.cells = 0, min.features = 0)

#Feature selection step: 
#Select highly variable genes with average expression set between 0.03 and 5, and z-score >= 5.
variable.genes <- FindVariableFeatures(seurat.data, selection.method = 'mvp', mean.cutoff = c(0.03, 5), dispersion.cutoff = c(0.5, Inf))
gene.names <- rownames(variable.genes)

#Scale data step:
#Shifts the mean expression of genes across cells to 0, and the variance of genes across cells to 1.
#Data is scaled according to highly variable genes so that these genes do not dominate during clustering.
scaled.data <- ScaleData(variable.genes, features = gene.names)

#Linear dimensional reduction step:
#PCA performed to filter out noise by extracting the first 50 PCs for guided clustering.
pca <- RunPCA(scaled.data, features = VariableFeatures(object = variable.genes))

#Cluster cells step:
clusters <- FindNeighbors(pca)
clusters <- FindClusters(clusters, resolution = 1)

