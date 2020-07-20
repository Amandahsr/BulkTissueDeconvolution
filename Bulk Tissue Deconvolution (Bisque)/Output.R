#This script details the output step for bulk tissue deconvolution.

#Extract cell-type proportions in each sample.
control.prop <- data.frame(control.est$bulk.props)

#Calculate mean of cell-type proportions across samples.
#Mean cell-type proportions are used as comparison with stained normal heart tissue to evaluate the accuracy of the Bisque algorithm.
control.avgprop <- data.frame(rowMeans(control.prop))
