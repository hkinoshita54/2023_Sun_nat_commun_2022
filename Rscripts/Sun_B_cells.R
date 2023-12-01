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
library(speckle)

# load data of B cells----
cts <- ReadMtx(mtx = "downloaded/B_cells/matrix.mtx.gz", features = "downloaded/B_cells/features.tsv.gz", cells = "downloaded/B_cells/barcodes.tsv.gz")
meta <- read.delim("downloaded/B_cells/meta_data.tsv", stringsAsFactors = F, header = T)
rownames(meta) <- meta[,1] ## rownames to cell id
meta <- meta[, -1] ## pick up necessary columns
seu <- CreateSeuratObject(counts = cts, meta.data = meta)
rm(cts, meta)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0)
table(seu$B_cluster)
table(seu$patient)

# cluster without integration----
# use counts slot as normalized data
seu <- NormalizeData(seu, verbose = FALSE)
seu[["RNA"]]$data <- seu[["RNA"]]$counts
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, cluster.name = "unintegrated_clusters", verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction.name = "umap.unintegrated", verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu, group.by = "batch") + NoAxes()
DimPlot(seu, group.by = "tissue") + NoAxes()
DimPlot(seu, group.by = "B_cluster") + NoAxes()

# cluster with scvi integration----
seu <- IntegrateLayers(
  object = seu, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)
seu <- FindNeighbors(seu, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 2, cluster.name = "scvi_clusters", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)
DimPlot(seu, reduction = "umap.scvi", label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu, reduction = "umap.scvi", group.by = "batch") + NoAxes()
DimPlot(seu, reduction = "umap.scvi", group.by = "tissue") + NoAxes()
DimPlot(seu, reduction = "umap.scvi", group.by = "B_cluster") + NoAxes()

saveRDS(seu, file = "RDSfiles/seu_B_downloaded.RDS")


