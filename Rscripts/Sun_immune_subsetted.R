# load packages----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
library(reticulate)
options(future.globals.maxSize = 1e10)
library(msigdbr)
library(DESeq2)
library(fgsea)

# load files----
seu <- readRDS("RDSfiles/seu_imm_subsetted.RDS")

# cluster
seu <- NormalizeData(seu, verbose = FALSE)
seu[["RNA"]]$data <- seu[["RNA"]]$counts
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- IntegrateLayers(
  object = seu, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)
seu <- FindNeighbors(seu, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.2, cluster.name = "scvi_clusters", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)
DimPlot(seu, reduction = "umap.scvi", label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu, reduction = "umap.scvi", group.by = "batch") + NoAxes()
DimPlot(seu, reduction = "umap.scvi", group.by = "tissue") + NoAxes()

FeaturePlot(seu,features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "IGHA1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "HBB", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# rename idents according to differentiation markers----
Idents(seu) <- "seurat_clusters"
seu <- RenameIdents(seu, 
                          `0` = "CD8Tcells", `1` = "CD4Tcells", `2` = "Bcells", `3` = "Treg", `4` = "NKcells", `5` = "CD4Tcells",`6` = "Plasma", `7` = "CD8Tcells", `8` = "Myeloids", `9` = "Mast", 
                          `10` = "CD8Tcells", `11` = "CD8Tcells", `12` = "CD8Tcells", `13` = "CD8Tcells", `14` = "NKcells", `15` = "Myeloids",`16` = "Erythrocytes", `17` = "Bcells", `18` = "Erythrocytes", `19` = "Myeloids"
                          )
DimPlot(seu, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
seu$cell_type <- Idents(seu)
seu$cell_type <- factor(x = seu$cell_type, levels = c("CD8Tcells", "CD4Tcells", "Treg", "NKcells", "Bcells", "Plasma", "Myeloids", "Mast", "Erythrocytes"))
saveRDS(seu, file = "RDSfiles/seu_imm_subsetted.RDS")

DimPlot(seu, reduction = "umap.scvi", group.by = "cluster") + NoAxes()
