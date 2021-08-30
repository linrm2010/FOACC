library(Seurat) 
library(dplyr) 
library(Matrix) 
library(patchwork) 
library(scales)

AthLeafSample.exp.filter <- Read10X(".")

AthLeafSample.exp.filter <- CreateSeuratObject(AthLeafSample.exp.filter, min.cells=3, min.features=200, project="AthLeafSample")

AthLeafSample.exp.filter[["percent.mt"]] <- PercentageFeatureSet(AthLeafSample.exp.filter, pattern = "^ATMG")
AthLeafSample.exp.filter[["percent.chl"]] <- PercentageFeatureSet(AthLeafSample.exp.filter, pattern = "^ATCG")

pdf("AthLeaf.step01.01.VlnPlot.pdf", height=6, width=8)
VlnPlot(AthLeafSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

AthLeafSample.exp.filter <- subset(AthLeafSample.exp.filter, subset = nFeature_RNA > 200 & nFeature_RNA < 8100 & percent.mt < 1 & percent.chl < 75)

pdf("AthLeaf.step01.01.VlnPlot.filter.pdf", height=6, width=8)
VlnPlot(AthLeafSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

AthLeafSample.exp.filter <- NormalizeData(AthLeafSample.exp.filter, normalization.method = "LogNormalize", scale.factor = 125000)

AthLeafSample.exp.filter <- FindVariableFeatures(AthLeafSample.exp.filter, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(AthLeafSample.exp.filter), 10)

pdf("AthLeaf.step01.02.VariableFeaturePlot.pdf", height=5, width=12)
VariableFeaturePlot(AthLeafSample.exp.filter) + LabelPoints(plot = VariableFeaturePlot(AthLeafSample.exp.filter), points = top10, repel = TRUE)
dev.off()

all.genes <- rownames(AthLeafSample.exp.filter)
AthLeafSample.exp.filter <- ScaleData(AthLeafSample.exp.filter, features = all.genes)

AthLeafSample.exp.filter <- RunPCA(AthLeafSample.exp.filter, npcs = 50, features = VariableFeatures(object = AthLeafSample.exp.filter))

pdf("AthLeaf.step01.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(AthLeafSample.exp.filter, dims = 1:2, reduction = "pca")
dev.off()

pdf("AthLeaf.step01.04.DimHeatmap.pdf", height=16, width=8)
DimHeatmap(AthLeafSample.exp.filter, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

AthLeafSample.exp.filter <- JackStraw(AthLeafSample.exp.filter, num.replicate = 100, dims = 50)
AthLeafSample.exp.filter <- ScoreJackStraw(AthLeafSample.exp.filter, dims = 1:30)

pdf("AthLeaf.step01.05.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(AthLeafSample.exp.filter, dims = 1:30)
dev.off()

pdf("AthLeaf.step01.06.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(AthLeafSample.exp.filter, ndims = 30)
dev.off()

AthLeafSample.exp.filter <- FindNeighbors(AthLeafSample.exp.filter, dims = 1:30)

AthLeafSample.exp.filter <- FindClusters(AthLeafSample.exp.filter, resolution = 0.5)

AthLeafSample.exp.filter.UMAP <- RunUMAP(AthLeafSample.exp.filter, dims = 1:30)

AthLeafSample.exp.filter.tSNE <- RunTSNE(AthLeafSample.exp.filter, dims = 1:30)

pdf("AthLeaf.step01.07.DimPlot.UMAP.pdf", height=3.5, width=5)
DimPlot(AthLeafSample.exp.filter.UMAP, reduction = "umap")
dev.off()

pdf("AthLeaf.step01.08.DimPlot.tSNE.pdf", height=3.5, width=5)
DimPlot(AthLeafSample.exp.filter.tSNE, reduction = "tsne")
dev.off()

color_list<-hue_pal()(13)
pdf("AthLeaf.step01.08.DimPlot.tSNE.select.pdf", height=3.2, width=4)
DimPlot(AthLeafSample.exp.filter.tSNE, reduction = "tsne", cols=color_list)
dev.off()

AthLeafSample.exp.filter.tSNE.markers <- FindAllMarkers(AthLeafSample.exp.filter.tSNE, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
AthLeafSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = myAUC)

top10 <- AthLeafSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)
pdf("AthLeaf.step01.09.DoHeatmap.pdf", height=8.5, width=8.5)
DoHeatmap(AthLeafSample.exp.filter.tSNE, features = top10$gene) + NoLegend()
dev.off()

saveRDS(AthLeafSample.exp.filter, file = "AthLeaf_step01.rds")

write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="0"), "Cells.step01.Ath01_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="1"), "Cells.step01.Ath02_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="2"), "Cells.step01.Ath03_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="3"), "Cells.step01.Ath04_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="4"), "Cells.step01.Ath05_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="5"), "Cells.step01.Ath06_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="6"), "Cells.step01.Ath07_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="7"), "Cells.step01.Ath08_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="8"), "Cells.step01.Ath09_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="9"), "Cells.step01.Ath10_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="10"), "Cells.step01.Ath11_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="11"), "Cells.step01.Ath12_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthLeafSample.exp.filter.tSNE, idents="12"), "Cells.step01.Ath13_ID.txt", sep="\t", quote = FALSE)

pdf("AthLeaf.step01.10.FeaturePlot.select.pdf", height=6.3, width=8)
FeaturePlot(AthLeafSample.exp.filter.tSNE, features = c("AT1G07640", "AT4G19840", "AT5G41920", "AT2G26250", "AT3G24140", "AT5G59870", "AT2G45190", "AT2G40100"))
dev.off()
