# UROPS: Bulk Tissue Cell-Type Deconvolution using scRNA-seq as Reference

This github details the code used in the UROPS project. 
The aim of this project is to evaluate the accuracy of bulk tissue deconvolution algorithms using scRNA-seq as reference.

The general project workflow is as follows:
Pre-processing of scRNA-seq data was first performed before guided clustering was carried out in Seurat. Bulk tissue deconvolution using the Bisque reference-based decomposition model was performed using 4 inputs: 1) unnormalised bulk RNA-seq, 2) unnormalised scRNA-seq, 3) cell-type proportions from guided clustering and 4) DEGs from guided clustering. Estimated cell-type proportions from deconvolution results were compared with stained normal heart tissues to evaluate the accuracy of cell-type proportions from Bisque.

Project overview:


