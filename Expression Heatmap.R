library(readxl)

#read dataset
dataset <- read_excel("~/Documents/Internships/UROPS/TPM.xlsx")
TPM <- dataset[,-1]

#appropriate data format. Added header "Genes" into excel TPM file to extract gene names for heatmap
tpm_matrix <- as.matrix(TPM)
rownames(tpm_matrix) <- dataset$Genes

#plot heatmap without clustering by rows or columns. Columns are samples, rows are genes.
heatmap(tpm_matrix, Colv = NA, Rowv = NA)