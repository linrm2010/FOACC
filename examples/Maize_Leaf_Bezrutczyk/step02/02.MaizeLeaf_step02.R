library(Seurat) 
library(dplyr) 
library(Matrix) 
library(patchwork) 
library(scales)

MaizeSample.exp.filter <- Read10X(".")

MaizeSample.exp.filter <- CreateSeuratObject(MaizeSample.exp.filter, min.cells=3, min.features=200, project="MaizeSample")

MaizeSample.exp.filter[["percent.mt"]] <- PercentageFeatureSet(MaizeSample.exp.filter, pattern = "^ZeamMp")
MaizeSample.exp.filter[["percent.chl"]] <- PercentageFeatureSet(MaizeSample.exp.filter, pattern = "^ZemaCp")

pdf("MaizeLeaf.step02.01.VlnPlot.pdf", height=6, width=8)
VlnPlot(MaizeSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

MaizeSample.exp.filter <- subset(MaizeSample.exp.filter, subset = nFeature_RNA > 200 & nFeature_RNA < 6700 & percent.mt < 2.5 & percent.chl < 20)

pdf("MaizeLeaf.step02.01.VlnPlot.filter.pdf", height=6, width=4)
VlnPlot(MaizeSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

MaizeSample.exp.filter <- NormalizeData(MaizeSample.exp.filter, normalization.method = "LogNormalize", scale.factor = 50000)

MaizeSample.exp.filter <- FindVariableFeatures(MaizeSample.exp.filter, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(MaizeSample.exp.filter), 10)

pdf("MaizeLeaf.step02.02.VariableFeaturePlot.pdf", height=5, width=12)
VariableFeaturePlot(MaizeSample.exp.filter) + LabelPoints(plot = VariableFeaturePlot(MaizeSample.exp.filter), points = top10, repel = TRUE)
dev.off()

all.genes <- rownames(MaizeSample.exp.filter)
MaizeSample.exp.filter <- ScaleData(MaizeSample.exp.filter, features = all.genes)

MaizeSample.exp.filter <- RunPCA(MaizeSample.exp.filter, npcs = 50, features = VariableFeatures(object = MaizeSample.exp.filter))

pdf("MaizeLeaf.step02.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(MaizeSample.exp.filter, dims = 1:2, reduction = "pca")
dev.off()

pdf("MaizeLeaf.step02.04.DimHeatmap.pdf", height=16, width=8)
DimHeatmap(MaizeSample.exp.filter, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

MaizeSample.exp.filter <- JackStraw(MaizeSample.exp.filter, num.replicate = 100, dims = 50)
MaizeSample.exp.filter <- ScoreJackStraw(MaizeSample.exp.filter, dims = 1:30)

pdf("MaizeLeaf.step02.05.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(MaizeSample.exp.filter, dims = 1:30)
dev.off()

pdf("MaizeLeaf.step02.06.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(MaizeSample.exp.filter, ndims = 30)
dev.off()

MaizeSample.exp.filter <- FindNeighbors(MaizeSample.exp.filter, dims = 1:16)

MaizeSample.exp.filter <- FindClusters(MaizeSample.exp.filter, resolution = 0.5)

MaizeSample.exp.filter.UMAP <- RunUMAP(MaizeSample.exp.filter, dims = 1:16)

MaizeSample.exp.filter.tSNE <- RunTSNE(MaizeSample.exp.filter, dims = 1:16)

pdf("MaizeLeaf.step02.07.DimPlot.UMAP.pdf", height=3.5, width=5)
DimPlot(MaizeSample.exp.filter.UMAP, reduction = "umap")
#DimPlot(MaizeSample.exp.filter.UMAP, reduction = "umap", label = TRUE)
dev.off()

pdf("MaizeLeaf.step02.08.DimPlot.tSNE.pdf", height=3.5, width=5)
DimPlot(MaizeSample.exp.filter.tSNE, reduction = "tsne")
#DimPlot(MaizeSample.exp.filter.tSNE, reduction = "tsne", label = TRUE)
dev.off()

color_list<-hue_pal()(10)
new_identity_TSNE<-MaizeSample.exp.filter.tSNE
new_identity_TSNE<-RenameIdents(new_identity_TSNE, "0"="2", "1"="0", "2"="3", "3"="1", "4"="5", "5"="4", "6"="6", "7"="7")
pdf("MaizeLeaf.step02.08.DimPlot.tSNE.select.pdf", height=3, width=3.7)
DimPlot(new_identity_TSNE, reduction = "tsne", cols=c(color_list[3], color_list[1], color_list[4], color_list[2], color_list[6], color_list[5], color_list[7], color_list[8]))
dev.off()
#      1 OriCluster      NewCluster	OriCluster	NewCluster
#     58 A08     B8	A08	B8
#    108 A07     B7	A07	B7
#    121 A05     B6	A05	B6
#    184 A06     B5	A06	B5
#    265 A02     B4	A02	B4
#    302 A04     B3	A04	B3
#    409 A01     B2	A01	B2
#    426 A03     B1	A03	B1

MaizeSample.exp.filter.tSNE.markers <- FindAllMarkers(MaizeSample.exp.filter.tSNE, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
MaizeSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = myAUC)

top10 <- MaizeSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)

pdf("MaizeLeaf.step02.09.DoHeatmap.pdf", height=8.5, width=8.5)
DoHeatmap(MaizeSample.exp.filter.tSNE, features = top10$gene) + NoLegend()
dev.off()

saveRDS(MaizeSample.exp.filter, file = "MaizeLeaf_step02.rds")

write.table(WhichCells(MaizeSample.exp.filter.tSNE, idents="0"), "Cells.step02.maize01.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(MaizeSample.exp.filter.tSNE, idents="1"), "Cells.step02.maize02.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(MaizeSample.exp.filter.tSNE, idents="2"), "Cells.step02.maize03.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(MaizeSample.exp.filter.tSNE, idents="3"), "Cells.step02.maize04.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(MaizeSample.exp.filter.tSNE, idents="4"), "Cells.step02.maize05.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(MaizeSample.exp.filter.tSNE, idents="5"), "Cells.step02.maize06.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(MaizeSample.exp.filter.tSNE, idents="6"), "Cells.step02.maize07.ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(MaizeSample.exp.filter.tSNE, idents="7"), "Cells.step02.maize08.ID.txt", sep="\t", quote = FALSE)
