library(data.table)
library(gplots)
library(ggplot2)

#Cleaning up data
data <- fread("~/Documents/Internships/UROPS/TPM.csv")
tpm <- data.frame(log10(data[,2:17]))
rownames(tpm) <- data$Genes

#Filter out highly variable genes
var.genes <- apply(tpm,1,sd)
var.names <- names(sort(var.genes, decreasing =TRUE))[1:500]  #select top 500 most variable genes
sub.genes <- tpm[var.names,]  #reduce dataset using variable gene names

#Pairwise Correlation between Samples
r.pairwise <- cor(sub.genes, method = "pearson")
cor.matrix <- as.matrix(1-r.pairwise)  #dissimilarity

#hierarchical clustering
hier.clus <- hclust(as.dist(1-r.pairwise), method="ward.D2")

#Expression Heatmap 
heatmap <- heatmap.2(cor.matrix, trace = "none", density.info = "none", Rowv=as.dendrogram(hier.clus))