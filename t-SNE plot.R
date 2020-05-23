library(readxl)
library(ggplot2)

#read dataset
dataset <- read_excel("~/Documents/Internships/UROPS/TPM.xlsx")
TPM <- dataset[,-1]

#appropriate data format
tpm_edited <- as.matrix(t(TPM))

#compute t-SNE
tsne <- Rtsne(tpm_edited, check_duplicates = FALSE, pca = TRUE, perplexity = 5, theta = 0.5, dims = 2)

#getting 2D matrix
tsne_matrix <- as.data.frame(tsne$Y) 

#plot t-SNE without clustering
tsne_plot <- ggplot(tsne_2D_matrix, aes(x=V1, y=V2))
tsne_plot + geom_point() + xlab("Dimension 1") + ylab("Dimension 2") + ggtitle("t-SNE plot") + theme_bw() + theme(panel.grid=element_blank())

#Not sure what perplexity in Rtsne function does, and also whether I should set pca argument in the function as TRUE or FALSE.