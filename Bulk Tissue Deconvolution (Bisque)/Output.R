#This script details the output step for bulk tissue deconvolution.
library(Bisque)
library(Biobase)

#Load deconvolution results.
control.est <- readRDS("path/control.est")

#Extract cell-type proportions from each sample.
control.prop <- data.frame(control.est$bulk.props)

#Calculate mean of cell-type proportions across samples. 
#Mean cell-type proportions are used as comparison with stained normal heart tissue to evaluate the accuracy of the Bisque algorithm.
control.avgprop <- data.frame(rowMeans(control.prop))
