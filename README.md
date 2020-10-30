# Heart Tissue Cell-Type Deconvolution using scRNA-seq as reference

This repository details the pipeline used to deconvolute bulk RNA-seq data using scRNA-seq as reference. 
The aim of this project is to evaluate the feasibility of applying reference-based bulk tissue deconvolution algorithms to heart tissues.

The general project workflow is as follows:
Pre-processing of scRNA-seq data was first performed before guided clustering was carried out in Seurat. Bulk tissue deconvolution using the Bisque reference-based decomposition model was performed using 4 inputs: 1) unnormalised bulk RNA-seq, 2) unnormalised scRNA-seq, 3) cell-type proportions from guided clustering and 4) DEGs from guided clustering. Estimated cell-type proportions from deconvolution results were compared with those from guided clustering to evaluate the accuracy of cell-type proportions from Bisque.

Project overview:
![](https://github.com/Amandahsr/UROPSBulkTissueDeconvolution/blob/master/Results/UROPS%20Project%20Overview.png)
