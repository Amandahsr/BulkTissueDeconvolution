library(readxl)
library(ggplot2)

#read dataset
dataset <- read_excel("~/Documents/Internships/UROPS/TPM.xlsx")
TPM <- dataset[,-1]

#transform into appropriate format
pca <- prcomp(t(TPM), scale = TRUE)
pca.data <- as.data.frame(pca.data)

#plot PCA using ggplot2
pca_plot <- ggplot(pca.data, aes(x=PC1, y=PC2))
pca_plot + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("PCA") + theme_bw() + theme(panel.grid=element_blank())

#%variance explained by each PC
percent_variance <- 100* pca$sdev^2/sum(pca$sdev^2)
round(percent_variance,2)

#plot %variance for each PC: scree plot
tpm.cov <- cov(pca.data)
tpm.eigen <- eigen(tpm.cov)
percent_variance <- tpm.eigen$values / sum(tpm.eigen$values)
tpm.screeplot <- qplot(c(1:16), percent_variance)
tpm.screeplot + geom_line() + xlab("PCs") + ylab("%Variance") + ggtitle("Scree Plot") + theme_bw() + theme(panel.grid=element_blank())

## How many PCs are you able to find?
#- 16 PCs, one for each column
## How do you determine the percentage variance explained by each PCs?
#-summary(pca), square the standard deviation then divide by the sum of all squared standard deviations
## Can you plot %variance explained by each PCs?
#- scree plot (?)