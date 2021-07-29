library(Seurat) 
library(dplyr) 
library(Matrix) 
library(patchwork) 
library(scales)

Cel.exp.filter <- Read10X(".")
Cel.exp.filter <- CreateSeuratObject(Cel.exp.filter, min.cells=3, min.features=200, project="10X_Cel")

pdf("Cel.step01.01.VlnPlot.pdf", height=6, width=8)
VlnPlot(Cel.exp.filter, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

Cel.exp.filter <- subset(Cel.exp.filter, subset = nFeature_RNA > 200 & nFeature_RNA < 2300)

pdf("Cel.step01.01.VlnPlot.filter.pdf", height=6, width=8)
VlnPlot(Cel.exp.filter, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

Cel.exp.filter <- NormalizeData(Cel.exp.filter, normalization.method = "LogNormalize", scale.factor = 6000)
Cel.exp.filter <- FindVariableFeatures(Cel.exp.filter, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(Cel.exp.filter), 10)
pdf("Cel.step01.02.VariableFeaturePlot.pdf", height=5, width=12)
VariableFeaturePlot(Cel.exp.filter) + LabelPoints(plot = VariableFeaturePlot(Cel.exp.filter), points = top10, repel = TRUE)
dev.off()

Cel.all.genes <- rownames(Cel.exp.filter)
Cel.exp.filter <- ScaleData(Cel.exp.filter, features = Cel.all.genes)

Cel.exp.filter <- RunPCA(Cel.exp.filter, npcs = 50, features = VariableFeatures(object = Cel.exp.filter))

pdf("Cel.step01.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(Cel.exp.filter, dims = 1:2, reduction = "pca")
dev.off()

pdf("Cel.step01.04.DimPlot.pdf", height=6, width=6)
DimPlot(Cel.exp.filter, reduction = "pca")
dev.off()

pdf("Cel.step01.05.DimHeatmap.pdf", height=16, width=8)
DimHeatmap(Cel.exp.filter, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

Cel.exp.filter <- JackStraw(Cel.exp.filter, num.replicate = 100, dims = 50)
Cel.exp.filter <- ScoreJackStraw(Cel.exp.filter, dims = 1:30)

pdf("Cel.step01.06.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(Cel.exp.filter, dims = 1:30)
dev.off()

pdf("Cel.step01.07.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(Cel.exp.filter, ndims = 30)
dev.off()

Cel.exp.filter <- FindNeighbors(Cel.exp.filter, dims = 1:30)

Cel.exp.filter <- FindClusters(Cel.exp.filter, resolution = 0.5)

Cel.exp.filter.UMAP <- RunUMAP(Cel.exp.filter, dims = 1:30)
Cel.exp.filter.tSNE <- RunTSNE(Cel.exp.filter, dims = 1:30)

pdf("Cel.step01.08.DimPlot.UMAP.pdf", height=3.5, width=5)
DimPlot(Cel.exp.filter.UMAP, reduction = "umap")
dev.off()

pdf("Cel.step01.09.DimPlot.tSNE.pdf", height=3.5, width=5)
DimPlot(Cel.exp.filter.tSNE, reduction = "tsne")
dev.off()

color_list<-hue_pal()(21)
pdf("Cel.step01.09.DimPlot.tSNE.select.pdf", height=5, width=6.5)
DimPlot(Cel.exp.filter.tSNE, reduction = "tsne", cols=color_list)
dev.off()

saveRDS(Cel.exp.filter, file = "Cel-Seurat4_step01.rds")

Cel.exp.filter.tSNE.markers <- FindAllMarkers(Cel.exp.filter.tSNE, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
Cel.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = myAUC)

top10 <- Cel.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)
pdf("Cel.step01.10.DoHeatmap.pdf", height=13, width=8)
DoHeatmap(Cel.exp.filter.tSNE, features = top10$gene) + NoLegend()
dev.off()

write.table(WhichCells(Cel.exp.filter.tSNE, idents="0"), "Cel.step01.Cell-00.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="1"), "Cel.step01.Cell-01.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="2"), "Cel.step01.Cell-02.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="3"), "Cel.step01.Cell-03.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="4"), "Cel.step01.Cell-04.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="5"), "Cel.step01.Cell-05.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="6"), "Cel.step01.Cell-06.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="7"), "Cel.step01.Cell-07.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="8"), "Cel.step01.Cell-08.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="9"), "Cel.step01.Cell-09.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="10"), "Cel.step01.Cell-10.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="11"), "Cel.step01.Cell-11.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="12"), "Cel.step01.Cell-12.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="13"), "Cel.step01.Cell-13.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="14"), "Cel.step01.Cell-14.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="15"), "Cel.step01.Cell-15.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="16"), "Cel.step01.Cell-16.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="17"), "Cel.step01.Cell-17.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="18"), "Cel.step01.Cell-18.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="19"), "Cel.step01.Cell-19.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="20"), "Cel.step01.Cell-20.ID.txt", sep="\t", quote = FALSE)
