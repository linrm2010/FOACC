library(Seurat) 
library(dplyr) 
library(Matrix) 
library(patchwork) 
library(scales)

RiceRootSample.exp.filter <- Read10X(".")

RiceRootSample.exp.filter <- CreateSeuratObject(RiceRootSample.exp.filter, min.cells=3, min.features=200, project="RiceRootSample")

pdf("RiceRoot.step02.01.VlnPlot.pdf", height=6, width=4)
VlnPlot(RiceRootSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

RiceRootSample.exp.filter <- subset(RiceRootSample.exp.filter, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)

pdf("RiceRoot.step02.01.VlnPlot.filter.pdf", height=6, width=4)
VlnPlot(RiceRootSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

RiceRootSample.exp.filter <- NormalizeData(RiceRootSample.exp.filter, normalization.method = "LogNormalize", scale.factor = 30000)

RiceRootSample.exp.filter <- FindVariableFeatures(RiceRootSample.exp.filter, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(RiceRootSample.exp.filter), 10)

pdf("RiceRoot.step02.02.VariableFeaturePlot.pdf", height=5, width=12)
VariableFeaturePlot(RiceRootSample.exp.filter) + LabelPoints(plot = VariableFeaturePlot(RiceRootSample.exp.filter), points = top10, repel = TRUE)
dev.off()

all.genes <- rownames(RiceRootSample.exp.filter)
RiceRootSample.exp.filter <- ScaleData(RiceRootSample.exp.filter, features = all.genes)

RiceRootSample.exp.filter <- RunPCA(RiceRootSample.exp.filter, npcs = 50, features = VariableFeatures(object = RiceRootSample.exp.filter))

pdf("RiceRoot.step02.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(RiceRootSample.exp.filter, dims = 1:2, reduction = "pca")
dev.off()

pdf("RiceRoot.step02.04.DimHeatmap.pdf", height=16, width=8)
DimHeatmap(RiceRootSample.exp.filter, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

RiceRootSample.exp.filter <- JackStraw(RiceRootSample.exp.filter, num.replicate = 100, dims = 50)
RiceRootSample.exp.filter <- ScoreJackStraw(RiceRootSample.exp.filter, dims = 1:30)

pdf("RiceRoot.step02.05.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(RiceRootSample.exp.filter, dims = 1:30)
dev.off()

pdf("RiceRoot.step02.06.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(RiceRootSample.exp.filter, ndims = 30)
dev.off()

RiceRootSample.exp.filter <- FindNeighbors(RiceRootSample.exp.filter, dims = 1:30)

RiceRootSample.exp.filter <- FindClusters(RiceRootSample.exp.filter, resolution = 0.5)

RiceRootSample.exp.filter.UMAP <- RunUMAP(RiceRootSample.exp.filter, dims = 1:30)

RiceRootSample.exp.filter.tSNE <- RunTSNE(RiceRootSample.exp.filter, dims = 1:30)

pdf("RiceRoot.step02.07.DimPlot.UMAP.pdf", height=3.5, width=5)
DimPlot(RiceRootSample.exp.filter.UMAP, reduction = "umap")
dev.off()

pdf("RiceRoot.step02.08.DimPlot.tSNE.pdf", height=3.5, width=5)
DimPlot(RiceRootSample.exp.filter.tSNE, reduction = "tsne")
dev.off()

color_list<-hue_pal()(14)
new_identity_TSNE<-RiceRootSample.exp.filter.tSNE
new_identity_TSNE<-RenameIdents(new_identity_TSNE, "0"="0+3", "1"="2", "2"="1", "3"="5", "4"="4", "5"="6", "6"="8", "7"="7+13", "8"="9", "9"="12", "10"="10", "11"="11")
pdf("RiceRoot.step02.08.DimPlot.tSNE.select.pdf", height=3.9, width=4.8)
DimPlot(new_identity_TSNE, reduction = "tsne", cols=c(color_list[1], color_list[3], color_list[2], color_list[6], color_list[5], color_list[7], color_list[9], color_list[8], color_list[10], color_list[13], color_list[11], color_list[12]))
dev.off()

#      1 OriCluster      NewCluster	OriCluster	NewCluster
#     17 A12     B12	A12	B12
#     26 A11     B11	A11	B11
#     58 A13     B10	A13	B10
#    140 A10     B09	A10	B09
#      8 A14     B08	A14	B08
#    144 A08     B08	A08	B08
#    154 A09     B07	A09	B07
#    228 A07     B06	A07	B06
#    287 A05     B05	A05	B05
#    299 A06     B04	A06	B04
#    398 A02     B03	A02	B03
#    591 A03     B02	A03	B02
#    565 A01     B01	A01	B01
#    138 A04     B01	A04	B01

RiceRootSample.exp.filter.tSNE.markers <- FindAllMarkers(RiceRootSample.exp.filter.tSNE, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
RiceRootSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = myAUC)

top10 <- RiceRootSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)
pdf("RiceRoot.step02.09.DoHeatmap.pdf", height=8.5, width=8.5)
DoHeatmap(RiceRootSample.exp.filter.tSNE, features = top10$gene) + NoLegend()
dev.off()

saveRDS(RiceRootSample.exp.filter, file = "RiceRoot_step02.rds")

write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="0"), "Cells.step02.rice01.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="1"), "Cells.step02.rice02.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="2"), "Cells.step02.rice03.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="3"), "Cells.step02.rice04.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="4"), "Cells.step02.rice05.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="5"), "Cells.step02.rice06.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="6"), "Cells.step02.rice07.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="7"), "Cells.step02.rice08.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="8"), "Cells.step02.rice09.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="9"), "Cells.step02.rice10.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="10"), "Cells.step02.rice11.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(RiceRootSample.exp.filter.tSNE, idents="11"), "Cells.step02.rice12.ID.txt", sep="\t", quote = FALSE)

pdf("RiceRoot.step02.10.FeaturePlot.select.pdf", height=5.5, width=9)
FeaturePlot(RiceRootSample.exp.filter.tSNE, features = c("LOC-Os03g25280", "LOC-Os03g12290", "LOC-Os04g46810", "LOC-Os08g03450", "LOC-Os03g61470", "LOC-Os03g37490", "LOC-Os07g07860", "LOC-Os06g38960", "LOC-Os01g73700", "LOC-Os01g73980", "LOC-Os10g42750"))
dev.off()
