setwd("E:/Spatial_Transcriptomics/Mouse Brain Coronal Section 2")

library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)

############## Data Loading ##################
brain <- Load10X_Spatial(data.dir = "E:/Spatial_Transcriptomics/Mouse Brain Coronal Section 2", filename = "CytAssist_FFPE_Mouse_Brain_Rep2_filtered_feature_bc_matrix.h5")

################# Data Processing #######################
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

############# Gene expression visualization ###########

SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

### Dimensionality reduction, clustering, and visualization ###
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

#visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with SpatialDimPlot().
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

### Interactive plotting ###
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
LinkedDimPlot(brain)

#### Identification of Spatially Variable Features #######
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
