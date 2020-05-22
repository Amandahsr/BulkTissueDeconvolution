library(readxl)
dataset <- read_excel("~/Documents/Internships/UROPS/TPM.xlsx")
TPM <- dataset[,-1]
pca <- prcomp(t(TPM), scale = TRUE)

library(ggplot2)
pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
pca.data <- as.data.frame(pca$x)
plot <- ggplot(pca.data, aes(x=PC1, y=PC2))
plot + geom_point()