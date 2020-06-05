library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)

#Read dataset
UMI <- fread("~/Documents/Internships/UROPS/GSE109816_normal_heart_umi_matrix.csv")
counts <- data.frame(UMI[,-1])
rownames(counts) <- UMI$V1
Sample.info <- fread("~/Documents/Internships/UROPS/GSE109816_normal_heart_cell_info.csv"
Labels <- cbind(Sample.info$ID, Sample.info$Type) 

#Merge UMI dataset with Tissue type labels
counts <- t(counts)
counts <- data.table(counts, keep.rownames = "ID")
merged.data <- merge(counts, Labels, by.x = "ID", by.y = "V1")

#Keep LV Samples, with counts >= 5 + >= 10 cells + top 30% genes
merged.data <- merged.data[which(merged.data$V2 != "N_LA_CM"),]
merged.data <- merged.data[which(merged.data$V2 != "N_LA_NCM"),]
reduced.data <- data.frame(apply(merged.data[,2:54751], 2, as.numeric))
rownames(reduced.data) <- merged.data$ID
reduced.data <- reduced.data[,colSums(reduced.data) >= 5] #genes with counts >= 5
reduced.data <- reduced.data[,colSums(reduced.data != 0) >= 10]  #genes with counts in >= 10 cells
reduced.data <- reduced.data[,rev(order(colSums(reduced.data)))]
reduced.data <- as.matrix(reduced.data[,1:8840])

#Normalise reduced data. Create Seurat object, sample names in columns, gene names in rows.
norm.data <- log(reduced.data+1,10)
seurat.data <- CreateSeuratObject(counts = t(norm.data), project = "umi.clustering", min.cells = 0, min.features = 0)

#Scale data
gene.names <- rownames(seurat.data)
scaled.data <- ScaleData(seurat.data, features = gene.names)

#PCA 
variable.data <- FindVariableFeatures(scaled.data)
variable.data <- RunPCA(variable.data, features = VariableFeatures(object = variable.data))

#Cluster cells and run t-SNE
tsne <- FindNeighbors(variable.data, dims = 1:10)
tsne <- FindClusters(tsne, resolution = 0.8)
tsne <- RunTSNE(tsne, perplexity = 5, check_duplicates = FALSE, reduction = "pca", seed.use = 42, tsne.method = "Rtsne", dim.embed = 2)
TSNEPlot(tsne, label = TRUE)

#Find DEGs in each cluster
markers <- FindAllMarkers(tsne, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(markers, "DEGs_UMI.csv") 
