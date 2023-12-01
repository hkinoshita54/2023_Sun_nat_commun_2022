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

# load data of all cells----
cts <- ReadMtx(mtx = "downloaded/all_cells/matrix.mtx.gz", features = "downloaded/all_cells/features.tsv.gz", cells = "downloaded/all_cells/barcodes.tsv.gz")
meta <- read.delim("downloaded/all_cells/meta_data.tsv", stringsAsFactors = F, header = T)
rownames(meta) <- meta[,1] ## rownames to cell id
meta <- meta[, -1] ## pick up necessary columns
seu <- CreateSeuratObject(counts = cts, meta.data = meta)
rm(cts, meta)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0)
table(seu$cluster)
table(seu$patient)

seu_epi <- subset(seu, cluster %in% c("Endocrine cells", "Epithelial cells"))
seu_str <- subset(seu, cluster %in% c("Endothelial cellls", "Fibroblasts", "Smooth muscle cells"))
seu_imm <- subset(seu, cluster %in% c("B cells", "B cells(Plasma cells)", "Erythrocytes", "Mast cells", "Myeloid cells", "T cells & NK cells"))
table(seu_epi$cluster)
table(seu_str$cluster)
table(seu_imm$cluster)

saveRDS(seu_epi, file = "RDSfiles/seu_epi_subsetted.RDS")
saveRDS(seu_str, file = "RDSfiles/seu_str_subsetted.RDS")
saveRDS(seu_imm, file = "RDSfiles/seu_imm_subsetted.RDS")

# load data of epithelial cells----
cts <- ReadMtx(mtx = "downloaded/epithelial_cells/matrix.mtx.gz", features = "downloaded/epithelial_cells/features.tsv.gz", cells = "downloaded/epithelial_cells/barcodes.tsv.gz")
meta <- read.delim("downloaded/epithelial_cells/meta_data.tsv", stringsAsFactors = F, header = T)
rownames(meta) <- meta[,1] ## rownames to cell id
meta <- meta[, -1] ## pick up necessary columns
seu <- CreateSeuratObject(counts = cts, meta.data = meta)
rm(cts, meta)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0)
table(seu$epi_cluster)
table(seu$patient)

saveRDS(seu, file = "RDSfiles/seu_epi_downloaded.RDS")

# compare seu_epi subsetted and seu from downloaded epithelial data----
seu_epi
seu
table(seu_epi$cluster)
table(seu$epi_cluster)
# endocrine cells are not included in downloaded epithelial data
# even after removing endocrine cells, subsetted seu_epi contains 4000 more cells

# cluster seu_epi without integration----
seu_epi <- NormalizeData(seu_epi, verbose = FALSE)
seu_epi <- FindVariableFeatures(seu_epi, verbose = FALSE)
seu_epi <- ScaleData(seu_epi, verbose = FALSE)
seu_epi <- RunPCA(seu_epi, npcs = 30, verbose = FALSE)
seu_epi <- FindNeighbors(seu_epi, dims = 1:30, verbose = FALSE)
seu_epi <- FindClusters(seu_epi, resolution = 1, cluster.name = "unintegrated_clusters", verbose = FALSE)
seu_epi <- RunUMAP(seu_epi, dims = 1:30, reduction.name = "umap.unintegrated", verbose = FALSE)
DimPlot(seu_epi, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu_epi, group.by = "batch") + NoAxes()
DimPlot(seu_epi, group.by = "tissue") + NoAxes()

FeaturePlot(seu_epi,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "ATP4B", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# too much immune contamination... maybe better to use downloaded epithelial data

# cluster seu from downloaded epithelial data without integration----
# use counts slot as normalized data
seu <- NormalizeData(seu, verbose = FALSE)
seu[["RNA"]]$data <- seu[["RNA"]]$counts
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, cluster.name = "unintegrated_clusters", verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction.name = "umap.unintegrated", verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu, group.by = "batch") + NoAxes()
DimPlot(seu, group.by = "tissue") + NoAxes()
DimPlot(seu, group.by = "epi_cluster") + NoAxes()

FeaturePlot(seu,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "ATP4B", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# no contamination of immune cells
# some of the samples make clusters by themselves as in the paper
# use this seurat object from downloaded epithelial data below

# Rank samples by ARID1A----
# use only tumor samples (remove paratumor samples)
seu_t <- subset(seu, tissue == "Tumor")
DimPlot(seu_t, group.by = "batch") + NoAxes()
DimPlot(seu_t, group.by = "epi_cluster") + NoAxes()
avg_t <- AverageExpression(seu_t, return.seurat = F, slot = "data", assays = "RNA", group.by = c("patient")) # use "data" to rank (not "count")
write.table(as.matrix(avg_t$RNA), "results/GeneExpression/avg_t_ver2.txt", sep="\t",col.names = T, row.names = T)  # check it in excel, rank samples by MUC6 levels
Idents(seu_t) <- "patient"
seu_t <- subset(seu_t, idents = "GC07", invert = T) # removed the median
ARID1A_H <- WhichCells(seu_t, ident = c("GC03",	"GC02",	"GC08",	"GC04"))
ARID1A_L <- WhichCells(seu_t, ident = c("GC06",	"GC05",	"GC09",	"GC10"))
seu_t <- SetIdent(seu_t, cells = ARID1A_H, value = "ARID1A_H")
seu_t$ARID1A_level <- Idents(seu_t)
seu_t <- SetIdent(seu_t, cells = ARID1A_L, value = "ARID1A_L")
seu_t$ARID1A_level <- Idents(seu_t)
seu_t$ARID1A_level <- factor(seu_t$ARID1A_level, levels = c("ARID1A_L", "ARID1A_H"))
DimPlot(seu_t, group.by = "ARID1A_level") + NoAxes()

# pseudobulk----
## need to be improved because DESeq2 should receive raw counts instead of processed counts
## maybe better to use other pipeline such as edgeR, limma etc.
bulk <- AggregateExpression(seu_t, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("patient", "ARID1A_level"))
tail(Cells(bulk))
bulk$patient <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$ARID1A_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$ARID1A_level <- factor(x = bulk$ARID1A_level, levels = c("L", "H"))
Idents(bulk) <- "ARID1A_level"
# bulk[["RNA"]]$counts <- as.matrix(bulk[["RNA"]]$counts) %>% round()
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)

write.table(as.matrix(de_markers),"./Results/DEG//DESeq2_Seurat_FindMarders_ARID1ALvsH_ver2.txt", sep ="\t", col.names = T,row.names = T)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val < 0.05, gene,"")), colour = "red", size = 3)

# do DEseq2 from matrix make rank for fgsea
cts <- as.matrix(bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
conditions <- factor(conditions, levels = c("L", "H"))
coldata <- cbind(colnames(cts), conditions) %>% as.data.frame()
coldata$conditions <- conditions
rownames(coldata) <- coldata$V1

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~conditions)
dds <- DESeq(dds)

vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conditions)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

res <- results(dds, contrast = c("conditions", "L", "H")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
write.table(as.matrix(res),"results/DEG//DESeq2_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, gene_name,"")), colour = "red", size = 3)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

# prepare gene sets
collections <- list()
collections$BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
collections$CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Homo sapiens", category = "H")
collections$C6 <- msigdbr(species = "Homo sapiens", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

# run fgsea
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_H_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)

plotEnrichment(collections$C6[["AKT_UP.V1_DN"]], ranks) + labs(title="C6_AKT_UP.V1_DN")
plotEnrichment(collections$C6[["AKT_UP.V1_UP"]], ranks) + labs(title="C6_AKT_UP.V1_UP")

# use DE genes from MUC6KO mouse RNAseq
ARID1AKO_UP_DN <- gmtPathways("gene_set/ARID1AKO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = ARID1AKO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(ARID1AKO_UP_DN[["ARID1AKO_UP"]], ranks) + labs(title="ARID1AKO_UP")
plotEnrichment(ARID1AKO_UP_DN[["ARID1AKO_DN"]], ranks) + labs(title="ARID1AKO_DN")

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_ARID1AKO_UP_DN_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)

# pseudobulk of fGSEA by ttest, rank by logFC*pvalue----
bulk <- AggregateExpression(seu_t, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("patient", "ARID1A_level"))
tail(Cells(bulk))
bulk$patient <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$ARID1A_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$ARID1A_level <- factor(x = bulk$ARID1A_level, levels = c("L", "H"))
Idents(bulk) <- "ARID1A_level"
# bulk[["RNA"]]$counts <- as.matrix(bulk[["RNA"]]$counts) %>% round()
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "t")
de_markers$gene <- rownames(de_markers)

write.table(as.matrix(de_markers),"./Results/DEG//ttest_Seurat_FindMarders_ARID1ALvsH_ver2.txt", sep ="\t", col.names = T,row.names = T)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val < 0.05, gene,"")), colour = "red", size = 3)

ranks <- de_markers$avg_log2FC * de_markers$p_val
names(ranks) <- rownames(de_markers)

# prepare gene sets
collections <- list()
collections$BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
collections$CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Homo sapiens", category = "H")
collections$C6 <- msigdbr(species = "Homo sapiens", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

# run fgsea
fgseaRes <- fgsea(pathways = collections$C6, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_C6_ARID1ALvsH_ttest_log2FC_pval.txt", sep ="\t", col.names = T,row.names = F)

plotEnrichment(collections$C6[["AKT_UP.V1_DN"]], ranks) + labs(title="C6_AKT_UP.V1_DN")
plotEnrichment(collections$C6[["AKT_UP.V1_UP"]], ranks) + labs(title="C6_AKT_UP.V1_UP")

# use DE genes from MUC6KO mouse RNAseq
ARID1AKO_UP_DN <- gmtPathways("gene_set/ARID1AKO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = ARID1AKO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(ARID1AKO_UP_DN[["ARID1AKO_UP"]], ranks) + labs(title="ARID1AKO_UP")
plotEnrichment(ARID1AKO_UP_DN[["ARID1AKO_DN"]], ranks) + labs(title="ARID1AKO_DN")

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_ARID1AKO_UP_DN_ARID1ALvsH_ttest_log2FC_pval.txt", sep ="\t", col.names = T,row.names = F)




















