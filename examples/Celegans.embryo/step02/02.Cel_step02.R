library(Seurat) 
library(dplyr) 
library(Matrix) 
library(patchwork) 
library(scales)

Cel.exp.filter <- Read10X(".")
Cel.exp.filter <- CreateSeuratObject(Cel.exp.filter, min.cells=3, min.features=200, project="10X_Cel")

pdf("Cel.step02.01.VlnPlot.pdf", height=6, width=8)
VlnPlot(Cel.exp.filter, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

Cel.exp.filter <- subset(Cel.exp.filter, subset = nFeature_RNA > 200 & nFeature_RNA < 2300)

pdf("Cel.step02.01.VlnPlot.filter.pdf", height=6, width=8)
VlnPlot(Cel.exp.filter, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

Cel.exp.filter <- NormalizeData(Cel.exp.filter, normalization.method = "LogNormalize", scale.factor = 6000)
Cel.exp.filter <- FindVariableFeatures(Cel.exp.filter, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(Cel.exp.filter), 10)
pdf("Cel.step02.02.VariableFeaturePlot.pdf", height=5, width=12)
VariableFeaturePlot(Cel.exp.filter) + LabelPoints(plot = VariableFeaturePlot(Cel.exp.filter), points = top10, repel = TRUE)
dev.off()

Cel.all.genes <- rownames(Cel.exp.filter)
Cel.exp.filter <- ScaleData(Cel.exp.filter, features = Cel.all.genes)

Cel.exp.filter <- RunPCA(Cel.exp.filter, npcs = 50, features = VariableFeatures(object = Cel.exp.filter))

pdf("Cel.step02.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(Cel.exp.filter, dims = 1:2, reduction = "pca")
dev.off()

pdf("Cel.step02.04.DimPlot.pdf", height=6, width=6)
DimPlot(Cel.exp.filter, reduction = "pca")
dev.off()

pdf("Cel.step02.05.DimHeatmap.pdf", height=16, width=10)
DimHeatmap(Cel.exp.filter, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

Cel.exp.filter <- JackStraw(Cel.exp.filter, num.replicate = 100, dims = 50)
Cel.exp.filter <- ScoreJackStraw(Cel.exp.filter, dims = 1:30)

pdf("Cel.step02.06.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(Cel.exp.filter, dims = 1:30)
dev.off()

pdf("Cel.step02.07.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(Cel.exp.filter, ndims = 30)
dev.off()

Cel.exp.filter <- FindNeighbors(Cel.exp.filter, dims = 1:20)

Cel.exp.filter <- FindClusters(Cel.exp.filter, resolution = 0.5)

Cel.exp.filter.UMAP <- RunUMAP(Cel.exp.filter, dims = 1:20)
Cel.exp.filter.tSNE <- RunTSNE(Cel.exp.filter, dims = 1:20)

pdf("Cel.step02.08.DimPlot.UMAP.pdf", height=5, width=6)
DimPlot(Cel.exp.filter.UMAP, reduction = "umap")
dev.off()

pdf("Cel.step02.09.DimPlot.tSNE.pdf", height=5, width=6)
DimPlot(Cel.exp.filter.tSNE, reduction = "tsne")
dev.off()

color_list<-hue_pal()(21) # 21, suggested by cluster numbers in step01
pdf("Cel.step02.09.DimPlot.tSNE.select.pdf", height=5, width=6)
DimPlot(Cel.exp.filter.tSNE, reduction = "tsne", cols=c(color_list[11], color_list[5], color_list[1], color_list[12], color_list[3], color_list[9], color_list[18], color_list[15], color_list[14], color_list[2], color_list[19], color_list[7], color_list[10], color_list[6], color_list[8], color_list[4], color_list[17], color_list[13], color_list[21]))
dev.off()
#      1 OriCluster      NewCluster	OriCluster	NewCluster
#     24 A21     B19	A21	B19
#     48 A13     B18	A13	B18
#     48 A17     B17	A17	B17
#     56 A04     B16	A04	B16
#     56 A08     B15	A08	B15
#     58 A06     B14	A06	B14
#     58 A10     B13	A10	B13
#     62 A07     B12	A07	B12
#     63 A19     B11	A19	B11
#     66 A02     B10	A02	B10
#     69 A14     B09	A14	B09
#     47 A15     B08	A15	B08
#     25 A16     B08	A16	B08
#     76 A18     B07	A18	B07
#     87 A09     B06	A09	B06
#     92 A03     B05	A03	B05
#    111 A12     B04	A12	B04
#    123 A01     B03	A01	B03
#    124 A05     B02	A05	B02
#     94 A11     B01	A11	B01
#     41 A20     B01	A20	B01

saveRDS(Cel.exp.filter, file = "Cel-Seurat4_step02.rds")

Cel.exp.filter.tSNE.markers <- FindAllMarkers(Cel.exp.filter.tSNE, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
Cel.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = myAUC)

top10 <- Cel.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)
pdf("Cel.step02.10.DoHeatmap.pdf", height=13, width=8)
DoHeatmap(Cel.exp.filter.tSNE, features = top10$gene) + NoLegend()
dev.off()

write.table(WhichCells(Cel.exp.filter.tSNE, idents="0"), "Cel.step02.Cell-00.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="1"), "Cel.step02.Cell-01.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="2"), "Cel.step02.Cell-02.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="3"), "Cel.step02.Cell-03.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="4"), "Cel.step02.Cell-04.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="5"), "Cel.step02.Cell-05.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="6"), "Cel.step02.Cell-06.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="7"), "Cel.step02.Cell-07.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="8"), "Cel.step02.Cell-08.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="9"), "Cel.step02.Cell-09.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="10"), "Cel.step02.Cell-10.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="11"), "Cel.step02.Cell-11.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="12"), "Cel.step02.Cell-12.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="13"), "Cel.step02.Cell-13.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="14"), "Cel.step02.Cell-14.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="15"), "Cel.step02.Cell-15.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="16"), "Cel.step02.Cell-16.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="17"), "Cel.step02.Cell-17.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(Cel.exp.filter.tSNE, idents="18"), "Cel.step02.Cell-18.ID.txt", sep="\t", quote = FALSE)
