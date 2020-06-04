library(data.table)
library(Rtsne)
library(ggplot2)

#read original datasets
UMI <- fread("~/Downloads/GSE109816_normal_heart_umi_matrix.csv")
counts <- UMI[,-1]   #remove gene names column
rownames(counts) <- UMI$V1
Sample.info <- fread("~/Downloads/GSE109816_normal_heart_cell_info.csv")
Labels <- cbind(Sample.info$ID, Sample.info$Type)      #retain only Sample ID and Tissue Type

#Merge UMI dataset with Tissue type labels
counts <- t(counts)   #Samples in rows, Genes in Columns
counts <- data.table(counts, keep.rownames = "ID")
merged.data <- merge(counts, Labels, by.x = "ID", by.y = "V1")

#Keep only LV Samples, with counts >= 5 & in >= 10 cells
merged.data <- merged.data[which(merged.data$V2.y != "N_LA_CM"),]
merged.data <- merged.data[which(merged.data$V2.y != "N_LA_NCM"),]
convert.data <- merged.data[, 2:54751]   #numeric data
convert.data <- apply(convert.data,2,as.numeric) 
convert.data <- data.table(t(convert.data))   #Samples in columns, Genes in rows
reduced.data <- convert.data[rowSums(convert.data) >= 5,]    #gene with counts >= 5 in at least 10 sample
reduced.data <- reduced.data[rowSums(reduced.data != 0)>=10,]    #gene with counts in at least 10 samples
reduced.data <- reduced.data[rev(order(rowSums(reduced.data))),]
reduced.data <- reduced.data[1:8840,]    #take top 30% highly expressed genes in each sample

#calulate tsne
normal.data <- log(reduced.data+1, 10)
set.seed(42)
tsne.data <- as.matrix(t(normal.data))  #Genes in columns, samples in rows
tsne <- Rtsne(tsne.data, check_duplicates = FALSE, pca = TRUE, perplexity = 5, theta = 0.5, dims = 2)
tsne.matrix <- as.data.frame(tsne$Y) 

#plot t-SNE
tsne_plot <- ggplot(tsne.matrix, aes(x=V1, y=V2, col=merged.data$V2.y))
tsne_plot + geom_point() + xlab("Dimension 1") + ylab("Dimension 2") + ggtitle("t-SNE Plot") + theme_bw()+ theme(panel.grid=element_blank()) + labs(colour = "Types")

