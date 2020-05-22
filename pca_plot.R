library(readxl)
library(ggplot2)

dataset <- read_excel("~/Documents/Internships/UROPS/TPM.xlsx")
TPM <- dataset[,-1]

pca <- prcomp(t(TPM), scale = TRUE)


pca.data <- as.data.frame(pca.data)
#data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
#pca.data <- as.data.frame(pca$x)

plot <- ggplot(pca.data, aes(x=PC1, y=PC2))
plot + geom_point() + theme_bw() + theme(panel.grid=element_blank())

## How many PCs are you able to find?
## How do you determine the percentage variance explained by each PCs?
## Can you plot %variance explained by each PCs?

