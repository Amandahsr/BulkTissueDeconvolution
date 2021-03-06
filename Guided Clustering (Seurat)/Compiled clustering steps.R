#This script details the full scRNA-seq guided clustering pipeline.
library(Seurat)
library(Matrix)

##1.Pre-processing steps:
#######################################################
#Read normalised 10X scRNA-seq dataset. The scRNA-seq data contains 287,269 cells from the 4 chambers of normal human heart samples (LV, LA, RV, RA) and 33,694 genes.
sc.counts <- Read10X("path/Single-Cell Data")

#Sample selection step: For the purpose of the project, retain only cells from the LV.
LV.cells <- sc.counts[,substr(colnames(sc.counts), 1, 2) == "LV"]

#Pre-process genes to ensure we only cluster genes with meaningful data. Retain genes that meet two criteria: 
#1) Has counts in at least 10 cells and 2) Belongs to the top 30% highly expressed genes.
LV.cells <- LV.cells[rowSums(LV.cells != 0) >= 10, ,drop = FALSE]

#Due to RAM limitations, only retain the top 30% highly-expressed genes.  
LV.cells <- LV.cells[rev(order(rowSums(LV.cells))),]
LV.cells <- LV.cells[1:7494,]

##2.Guided clustering steps:
#######################################################
#Initialise UMI counts as Seurat object. min.cells and min.features are set to 0 because pre-processing was already performed.
seurat.data <- CreateSeuratObject(LV.cells, project = 'Seurat Guided Clustering', min.cells = 0, min.features = 0)

#Select highly variable genes with average expression set between 0.03 and 5, and z-score >= 5. This highlights meaningful biological signals during clustering.
variable.genes <- FindVariableFeatures(seurat.data, selection.method = 'mvp', mean.cutoff = c(0.03, 5), dispersion.cutoff = c(0.5, Inf))
gene.names <- rownames(variable.genes)

#Shift the mean expression of genes across cells to 0, and the variance of genes across cells to 1. Data is scaled according to highly variable genes so that these genes do not dominate during clustering.
scaled.data <- ScaleData(variable.genes, features = gene.names)

#Linear dimensional reduction step: PCA is performed to filter out noise by extracting the first 50 PCs for guided clustering.
pca <- RunPCA(scaled.data, features = VariableFeatures(object = variable.genes))

#Cluster cells by creating a KNN graph and optimizing it.
clusters <- FindNeighbors(pca)
clusters <- FindClusters(clusters, resolution = 1)

##3.Output steps:
#######################################################
#UMAP: Visualise cell clusters using only first 50 PCs.
umap <- RunUMAP(clusters, dims = 1:50, n.neighbors = 25, min.dist = 1.5, spread = 3.5)
DimPlot(umap, reduction = "umap", label = TRUE)

#DEG analysis: Find genes that are preferentially expressed in each cluster with p-value <= 0.05.
clusters.markers <- FindAllMarkers(umap, only.pos=TRUE, return.thresh = 0.05)

#Feature plots: Visualise cell-type markers to identify cell types. A subset of markers used to identify cell-types are shown below.
Cardiomyocyte.markers <- FeaturePlot(umap, features = c('MYH7', 'TTN', 'TNNT2', 'RYR2'))
Macrophage.markers <- FeaturePlot(umap, features = c('CD163', 'MRC1', 'COLEC12', 'MARCH1', 'SLC11A1', 'RBPJ', 'F13A1'))
SmoothMuscleCell.markers <- FeaturePlot(umap, features = c('MYH11', 'LMOD1'))
Fibroblast.markers <- FeaturePlot(umap, features = c('POSTN', 'IGF1', 'ADAMTS4', 'DCN', 'ELN', 'NOX4', 'VCAN'))
Endothelial.markers <- FeaturePlot(umap, features = c('PECAM1', 'PROX1', 'VWF'))
Pericyte.markers <- FeaturePlot(umap, features = c('PDGFRB', 'STEAP4'))
TCell.markers <- FeaturePlot(umap, features = c('PTPRC', 'SKAP1'))
BCell.markers <- FeaturePlot(umap, features = c('BANK1'))
NeuronalCell.markers <- FeaturePlot(umap, features = c('NRXN1', 'NRXN3', 'NCAM2'))

#Relabel UMAP using cell coordinates after performing pathway analysis. Cluster 24 was found to contain two distinct cell-types, and is relabelled as two distinct clusters 24a and 24b.
umap.coordinates <- FetchData(umap, vars = c("ident", "UMAP_1", "UMAP_2"))
umap.relabelled <- umap.coordinates[(umap.coordinates$ident == '24'),]
TCell.cluster <- umap.relabelled[umap.relabelled$UMAP_1 > 0,]
NeuronalCell.cluster <- umap.relabelled[umap.relabelled$UMAP_1 < 0,]
Idents(umap, cells = rownames(TCell.cluster)) <- '24a'
Idents(umap, cells = rownames(NeuronalCell.cluster)) <- '24b'
DimPlot(umap, reduction = "umap", label = TRUE)

#Relabel clusters by cell types
CM.cells <- umap.coordinates[(umap.coordinates$ident == '0')|(umap.coordinates$ident == '2')|
                               (umap.coordinates$ident == '4')|(umap.coordinates$ident == '8')|
                               (umap.coordinates$ident == '9')|(umap.coordinates$ident == '10')|
                               (umap.coordinates$ident == '12')|(umap.coordinates$ident == '16')|
                               (umap.coordinates$ident == '18')|(umap.coordinates$ident == '20')|
                               (umap.coordinates$ident == '23')|(umap.coordinates$ident == '27'),]
Idents(umap, cells = rownames(CM.cells)) <- 'Cardiomyocytes'
FC.cells <- umap.coordinates[(umap.coordinates$ident == '1')|(umap.coordinates$ident == '6')|
                               (umap.coordinates$ident == '7')|(umap.coordinates$ident == '13')|
                               (umap.coordinates$ident == '15')|(umap.coordinates$ident == '21')|
                               (umap.coordinates$ident == '28'),]
Idents(umap, cells = rownames(FC.cells)) <- 'Fibroblasts'
EC.cells <- umap.coordinates[(umap.coordinates$ident == '5')|(umap.coordinates$ident == '14')|
                               (umap.coordinates$ident == '19')|(umap.coordinates$ident == '25')|
                               (umap.coordinates$ident == '29'),]
Idents(umap, cells = rownames(EC.cells)) <- 'Endothelial Cells'
SMC.cells <- umap.coordinates[(umap.coordinates$ident == '3')|(umap.coordinates$ident == '30'),]
Idents(umap, cells = rownames(SMC.cells)) <- 'Smooth Muscle Cells'
MC.cells <- umap.coordinates[(umap.coordinates$ident == '11')|(umap.coordinates$ident == '22'),]
Idents(umap, cells = rownames(MC.cells)) <- 'Macrophages'
PC.cells <- umap.coordinates[(umap.coordinates$ident == '17'),]
Idents(umap, cells = rownames(PC.cells)) <- 'Pericytes'
Idents(umap, cells = rownames(TCell.cluster)) <- 'T-Cells'
Idents(umap, cells = rownames(NeuronalCell.cluster)) <- 'Neuronal Cells'
BC.cells <- umap.coordinates[(umap.coordinates$ident == '26'),]
Idents(umap, cells = rownames(BC.cells)) <- 'B-Cells'
DimPlot(umap, reduction = "umap", label = TRUE)

#Perform second iteration of DEG analysis after relabelling to confirm cell types.
clusters2.markers <- FindAllMarkers(umap, only.pos=TRUE, return.thresh = 0.05)

#Extract cell-type proportions as input for bulk-tissue deconvolution.
celltype.proportions <- prop.table(table(Idents(umap)))

