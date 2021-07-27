library(Seurat) 
library(dplyr) 
library(Matrix) 
library(patchwork) 
library(scales)

AthRootSample.exp.filter <- Read10X(".")

AthRootSample.exp.filter <- CreateSeuratObject(AthRootSample.exp.filter, min.cells=3, min.features=200, project="AthRootSample")

AthRootSample.exp.filter[["percent.mt"]] <- PercentageFeatureSet(AthRootSample.exp.filter, pattern = "^ATMG")
AthRootSample.exp.filter[["percent.chl"]] <- PercentageFeatureSet(AthRootSample.exp.filter, pattern = "^ATCG")

pdf("AthRoot.step01.01.VlnPlot.pdf", height=6, width=8)
VlnPlot(AthRootSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

AthRootSample.exp.filter <- subset(AthRootSample.exp.filter, subset = nFeature_RNA > 200 & nFeature_RNA < 9600 & percent.mt < 0.07 & percent.chl < 0.1)

pdf("AthRoot.step01.01.VlnPlot.filter.pdf", height=6, width=8)
VlnPlot(AthRootSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

AthRootSample.exp.filter <- NormalizeData(AthRootSample.exp.filter, normalization.method = "LogNormalize", scale.factor = 150000)

AthRootSample.exp.filter <- FindVariableFeatures(AthRootSample.exp.filter, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(AthRootSample.exp.filter), 10)

pdf("AthRoot.step01.02.VariableFeaturePlot.pdf", height=5, width=12)
VariableFeaturePlot(AthRootSample.exp.filter) + LabelPoints(plot = VariableFeaturePlot(AthRootSample.exp.filter), points = top10, repel = TRUE)
dev.off()

all.genes <- rownames(AthRootSample.exp.filter)
AthRootSample.exp.filter <- ScaleData(AthRootSample.exp.filter, features = all.genes)

AthRootSample.exp.filter <- RunPCA(AthRootSample.exp.filter, npcs = 50, features = VariableFeatures(object = AthRootSample.exp.filter))

pdf("AthRoot.step01.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(AthRootSample.exp.filter, dims = 1:2, reduction = "pca")
dev.off()

pdf("AthRoot.step01.04.DimPlot.pdf", height=6, width=6)
DimPlot(AthRootSample.exp.filter, reduction = "pca")
dev.off()

pdf("AthRoot.step01.05.DimHeatmap.pdf", height=16, width=8)
DimHeatmap(AthRootSample.exp.filter, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

AthRootSample.exp.filter <- JackStraw(AthRootSample.exp.filter, num.replicate = 100, dims = 50)
AthRootSample.exp.filter <- ScoreJackStraw(AthRootSample.exp.filter, dims = 1:30)

pdf("AthRoot.step01.06.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(AthRootSample.exp.filter, dims = 1:30)
dev.off()

pdf("AthRoot.step01.07.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(AthRootSample.exp.filter, ndims = 30)
dev.off()

AthRootSample.exp.filter <- FindNeighbors(AthRootSample.exp.filter, dims = 1:30)

AthRootSample.exp.filter <- FindClusters(AthRootSample.exp.filter, resolution = 0.5)

AthRootSample.exp.filter.UMAP <- RunUMAP(AthRootSample.exp.filter, dims = 1:30)
AthRootSample.exp.filter.tSNE <- RunTSNE(AthRootSample.exp.filter, dims = 1:30)

pdf("AthRoot.step01.08.DimPlot.UMAP.pdf", height=3.5, width=5)
DimPlot(AthRootSample.exp.filter.UMAP, reduction = "umap")
dev.off()

pdf("AthRoot.step01.09.DimPlot.tSNE.pdf", height=3.5, width=5)
DimPlot(AthRootSample.exp.filter.tSNE, reduction = "tsne")
dev.off()

AthRootSample.exp.filter.tSNE.markers <- FindAllMarkers(AthRootSample.exp.filter.tSNE, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
AthRootSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = myAUC)

top10 <- AthRootSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)

pdf("AthRoot.step01.10.DoHeatmap.pdf", height=8, width=8)
DoHeatmap(AthRootSample.exp.filter.tSNE, features = top10$gene) + NoLegend()
dev.off()

saveRDS(AthRootSample.exp.filter, file = "AthRoot_step01.rds")

write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="0"), "AthRoot.step01.Ath01_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="1"), "AthRoot.step01.Ath02_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="2"), "AthRoot.step01.Ath03_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="3"), "AthRoot.step01.Ath04_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="4"), "AthRoot.step01.Ath05_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="5"), "AthRoot.step01.Ath06_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="6"), "AthRoot.step01.Ath07_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="7"), "AthRoot.step01.Ath08_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="8"), "AthRoot.step01.Ath09_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="9"), "AthRoot.step01.Ath10_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="10"), "AthRoot.step01.Ath11_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="11"), "AthRoot.step01.Ath12_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="12"), "AthRoot.step01.Ath13_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="13"), "AthRoot.step01.Ath14_ID.txt", sep="\t", quote = FALSE)
