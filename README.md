# Spatial Transcriptomics Analysis with Seurat and Shiny

## Overview
This project provides an analysis pipeline for spatial transcriptomics data using the **Seurat** package in R. It processes and visualizes gene expression data from a **Mouse Brain Coronal Section** dataset.

## Prerequisites
Ensure you have the following R packages installed before running the script:
- **shiny**
- **Seurat**
- **ggplot2**
- **patchwork**

You can install these packages using the following commands:
```r
install.packages("shiny")
install.packages("ggplot2")
install.packages("patchwork")
if (!requireNamespace("Seurat", quietly = TRUE)) {
    install.packages("Seurat")
}
```

## Data Preparation
The analysis is based on a **10X Genomics Spatial Transcriptomics dataset** stored in an HDF5 file format. Ensure the dataset is placed in the correct directory:
```
E:/Spatial_Transcriptomics/Mouse Brain Coronal Section 2
```
The dataset should contain the file:
```
CytAssist_FFPE_Mouse_Brain_Rep2_filtered_feature_bc_matrix.h5
```

## Workflow

### 1. Load Data
The dataset is loaded using the `Load10X_Spatial()` function from Seurat:
```r
brain <- Load10X_Spatial(data.dir = "E:/Spatial_Transcriptomics/Mouse Brain Coronal Section 2",
                         filename = "CytAssist_FFPE_Mouse_Brain_Rep2_filtered_feature_bc_matrix.h5")
```

### 2. Data Processing
- Visualization of spatial feature counts:
```r
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
```
- Normalization and transformation using SCTransform:
```r
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
```

### 3. Gene Expression Visualization
- Visualizing specific genes:
```r
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
```

### 4. Clustering and Dimensionality Reduction
- Perform PCA, find clusters, and run UMAP:
```r
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
```
- Visualization of clusters:
```r
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
```

### 5. Interactive Visualization
- Enable interactive plotting for gene expression:
```r
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
LinkedDimPlot(brain)
```

### 6. Identification of Spatially Variable Features
- Find differentially expressed markers between clusters:
```r
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
```

## Usage
Run the R script in an environment that supports Shiny and Seurat, such as RStudio. Ensure the dataset is in the correct directory before execution.

## Acknowledgments
This analysis leverages **Seurat** for processing spatial transcriptomics data. Special thanks to the developers of **10X Genomics** and the **Seurat** package for enabling spatial transcriptomics research.

---

For questions or issues, please reach out via GitHub or relevant scientific forums.

