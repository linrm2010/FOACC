library(Seurat)
library(dplyr)
library(patchwork)

pbmc.data <- Read10X(data.dir = ".")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pdf("pbmc.step01.01.VlnPlot.pdf", height=3.5, width=8)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("pbmc.step01.01.FeatureScatter.pdf", height=3.5, width=8)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)

pdf("pbmc.step01.02.VariableFeaturePlot.pdf", height=3, width=10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pdf("pbmc.step01.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
dev.off()

pdf("pbmc.step01.04.DimPlot.pdf", height=6, width=6)
DimPlot(pbmc, reduction = "pca")
dev.off()

pdf("pbmc.step01.05.DimHeatmap.pdf", height=16, width=8)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

pdf("pbmc.step01.06.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(pbmc, dims = 1:15)
dev.off()

pdf("pbmc.step01.07.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(pbmc)
dev.off()

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc.UMAP <- RunUMAP(pbmc, dims = 1:10)
pbmc.tSNE <- RunTSNE(pbmc, dims = 1:10)

pdf("pbmc.step01.08.DimPlot.UMAP.pdf", height=4, width=5)
DimPlot(pbmc.UMAP, reduction = "umap")
dev.off()

pdf("pbmc.step01.09.DimPlot.tSNE.pdf", height=4, width=5)
DimPlot(pbmc.tSNE, reduction = "tsne")
dev.off()

library(scales)
color_list<-hue_pal()(9)
pdf("pbmc.step01.09.DimPlot.tSNE.select.pdf", height=3.2, width=4)
DimPlot(pbmc.tSNE, reduction = "tsne", cols=color_list)
dev.off()

write.table(WhichCells(pbmc.UMAP, idents="0"), "pbmc_UMAP_cell.0.txt", sep="\t", quote = FALSE)
write.table(WhichCells(pbmc.UMAP, idents="1"), "pbmc_UMAP_cell.1.txt", sep="\t", quote = FALSE)
write.table(WhichCells(pbmc.UMAP, idents="2"), "pbmc_UMAP_cell.2.txt", sep="\t", quote = FALSE)
write.table(WhichCells(pbmc.UMAP, idents="3"), "pbmc_UMAP_cell.3.txt", sep="\t", quote = FALSE)
write.table(WhichCells(pbmc.UMAP, idents="4"), "pbmc_UMAP_cell.4.txt", sep="\t", quote = FALSE)
write.table(WhichCells(pbmc.UMAP, idents="5"), "pbmc_UMAP_cell.5.txt", sep="\t", quote = FALSE)
write.table(WhichCells(pbmc.UMAP, idents="6"), "pbmc_UMAP_cell.6.txt", sep="\t", quote = FALSE)
write.table(WhichCells(pbmc.UMAP, idents="7"), "pbmc_UMAP_cell.7.txt", sep="\t", quote = FALSE)
write.table(WhichCells(pbmc.UMAP, idents="8"), "pbmc_UMAP_cell.8.txt", sep="\t", quote = FALSE)

pbmc.markers <- FindAllMarkers(pbmc.UMAP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

pdf("pbmc.step01.10.VlnPlot.pdf", height=4, width=8)
VlnPlot(pbmc.UMAP, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
dev.off()

pdf("pbmc.step01.11.FeaturePlot.pdf", height=7, width=8)
FeaturePlot(pbmc.UMAP, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
dev.off()

# MS4A1 <-> ENSG00000156738, GNLY <-> ENSG00000115523, CD3E <-> ENSG00000198851
# CD14 <-> ENSG00000170458, FCER1A <-> ENSG00000179639, FCGR3A <-> ENSG00000203747
# LYZ <-> ENSG00000090382, PPBP <-> ENSG00000163736, CD8A <-> ENSG00000153563
# FeaturePlot(pbmc.tSNE, features = c("ENSG00000156738", "ENSG00000115523", "ENSG00000198851", "ENSG00000170458", "ENSG00000179639", "ENSG00000203747", "ENSG00000090382", "ENSG00000163736", "ENSG00000153563"))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("pbmc.step01.12.DoHeatmap.pdf", height=8, width=8)
DoHeatmap(pbmc.UMAP, features = top10$gene) + NoLegend()
dev.off()

# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc.UMAP)
# pbmc <- RenameIdents(pbmc.UMAP, new.cluster.ids)
# DimPlot(pbmc.UMAP, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "pbmc3k_step01.rds")
