library(readxl)
library(ggplot2)

<<<<<<< HEAD
#read dataset
dataset <- read_excel("~/Documents/Internships/UROPS/TPM.xlsx")
TPM <- dataset[,-1]

#transform into appropriate format
=======
dataset <- read_excel("~/Documents/Internships/UROPS/TPM.xlsx")
TPM <- dataset[,-1]

>>>>>>> 5fbe7555bfb2c95ae0de34e585272e68e393f1ed
pca <- prcomp(t(TPM), scale = TRUE)
pca.data <- as.data.frame(pca.data)

<<<<<<< HEAD
#plot PCA using ggplot2
pca_plot <- ggplot(pca.data, aes(x=PC1, y=PC2))
pca_plot + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("PCA") + theme_bw() + theme(panel.grid=element_blank())

#plot %variance for each PC: scree plot
tpm.cov <- cov(pca.data)
tpm.eigen <- eigen(tpm.cov)
percent_variance <- tpm.eigen$values / sum(tpm.eigen$values)
tpm.screeplot <- qplot(c(1:16), percent_variance)
tpm.screeplot + geom_line() + xlab("PCs") + ylab("%Variance") + ggtitle("Scree Plot") + theme_bw() + theme(panel.grid=element_blank())

## How many PCs are you able to find? - 16 PCs, one for each column
## How do you determine the percentage variance explained by each PCs? -summary(pca) gives this under 'proportion of variance'
## Can you plot %variance explained by each PCs? - scree plot
=======

pca.data <- as.data.frame(pca.data)
#data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
#pca.data <- as.data.frame(pca$x)

plot <- ggplot(pca.data, aes(x=PC1, y=PC2))
plot + geom_point() + theme_bw() + theme(panel.grid=element_blank())

## How many PCs are you able to find?
## How do you determine the percentage variance explained by each PCs?
## Can you plot %variance explained by each PCs?

>>>>>>> 5fbe7555bfb2c95ae0de34e585272e68e393f1ed
