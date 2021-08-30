library(Seurat) 
library(dplyr) 
library(Matrix) 
library(patchwork) 
library(scales)

RiceRootSample.exp.filter <- Read10X(".")

RiceRootSample.exp.filter <- CreateSeuratObject(RiceRootSample.exp.filter, min.cells=3, min.features=200, project="RiceRootSample")

pdf("RiceRoot.step01.01.VlnPlot.pdf", height=6, width=8)
VlnPlot(RiceRootSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

RiceRootSample.exp.filter <- subset(RiceRootSample.exp.filter, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)

pdf("RiceRoot.step01.01.VlnPlot.filter.pdf", height=6, width=4)
VlnPlot(RiceRootSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

RiceRootSample.exp.filter <- NormalizeData(RiceRootSample.exp.filter, normalization.method = "LogNormalize", scale.factor = 40000)

RiceRootSample.exp.filter <- FindVariableFeatures(RiceRootSample.exp.filter, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(RiceRootSample.exp.filter), 10)

pdf("RiceRoot.step01.02.VariableFeaturePlot.pdf", height=5, width=12)
VariableFeaturePlot(RiceRootSample.exp.filter) + LabelPoints(plot = VariableFeaturePlot(RiceRootSample.exp.filter), points = top10, repel = TRUE)
dev.off()

all.genes <- rownames(RiceRootSample.exp.filter)
RiceRootSample.exp.filter <- ScaleData(RiceRootSample.exp.filter, features = all.genes)

RiceRootSample.exp.filter <- RunPCA(RiceRootSample.exp.filter, npcs = 50, features = VariableFeatures(object = RiceRootSample.exp.filter))

pdf("RiceRoot.step01.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(RiceRootSample.exp.filter, dims = 1:2, reduction = "pca")
dev.off()

pdf("RiceRoot.step01.04.DimHeatmap.pdf", height=16, width=8)
DimHeatmap(RiceRootSample.exp.filter, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

RiceRootSample.exp.filter <- JackStraw(RiceRootSample.exp.filter, num.replicate = 100, dims = 50)
RiceRootSample.exp.filter <- ScoreJackStraw(RiceRootSample.exp.filter, dims = 1:30)

pdf("RiceRoot.step01.05.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(RiceRootSample.exp.filter, dims = 1:30)
dev.off()

pdf("RiceRoot.step01.06.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(RiceRootSample.exp.filter, ndims = 30)
dev.off()

RiceRootSample.exp.filter <- FindNeighbors(RiceRootSample.exp.filter, dims = 1:30)

RiceRootSample.exp.filter <- FindClusters(RiceRootSample.exp.filter, resolution = 0.5)

RiceRootSample.exp.filter.UMAP <- RunUMAP(RiceRootSample.exp.filter, dims = 1:30)

RiceRootSample.exp.filter.tSNE <- RunTSNE(RiceRootSample.exp.filter, dims = 1:30)

pdf("RiceRoot.step01.07.DimPlot.UMAP.pdf", height=3.5, width=5)
DimPlot(RiceRootSample.exp.filter.UMAP, reduction = "umap")
dev.off()

pdf("RiceRoot.step01.08.DimPlot.tSNE.pdf", height=3.5, width=5)
DimPlot(RiceRootSample.exp.filter.tSNE, reduction = "tsne")
dev.off()

color_list<-hue_pal()(14)
pdf("RiceRoot.step01.08.DimPlot.tSNE.select.pdf", height=3.9, width=4.8)
DimPlot(RiceRootSample.exp.filter.tSNE, reduction = "tsne", cols=color_list)
dev.off()

RiceRootSample.exp.filter.tSNE.markers <- FindAllMarkers(RiceRootSample.exp.filter.tSNE, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
RiceRootSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = myAUC)

top10 <- RiceRootSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)
pdf("RiceRoot.step01.09.DoHeatmap.pdf", height=8.5, width=8.5)
DoHeatmap(RiceRootSample.exp.filter.tSNE, features = top10$gene) + NoLegend()
dev.off()

saveRDS(RiceRootSample.exp.filter, file = "RiceRoot_step01.rds")

write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="0"), "Cells.step01.rice01.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="1"), "Cells.step01.rice02.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="2"), "Cells.step01.rice03.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="3"), "Cells.step01.rice04.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="4"), "Cells.step01.rice05.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="5"), "Cells.step01.rice06.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="6"), "Cells.step01.rice07.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="7"), "Cells.step01.rice08.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="8"), "Cells.step01.rice09.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="9"), "Cells.step01.rice10.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="10"), "Cells.step01.rice11.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="11"), "Cells.step01.rice12.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="12"), "Cells.step01.rice13.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="13"), "Cells.step01.rice14.ID.txt", sep="\t", quote = FALSE)

pdf("RiceRoot.step01.10.FeaturePlot.select.pdf", height=5.5, width=9)
FeaturePlot(RiceRootSample.exp.filter.tSNE, features = c("LOC-Os03g25280", "LOC-Os03g12290", "LOC-Os04g46810", "LOC-Os08g03450", "LOC-Os03g61470", "LOC-Os03g37490", "LOC-Os07g07860", "LOC-Os06g38960", "LOC-Os01g73700", "LOC-Os01g73980", "LOC-Os10g42750"))
dev.off()
