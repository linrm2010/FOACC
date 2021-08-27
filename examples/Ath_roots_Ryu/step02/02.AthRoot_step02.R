library(Seurat) 
library(dplyr) 
library(Matrix) 
library(patchwork) 
library(scales)

AthRootSample.exp.filter <- Read10X(".")

AthRootSample.exp.filter <- CreateSeuratObject(AthRootSample.exp.filter, min.cells=3, min.features=200, project="AthRootSample")

AthRootSample.exp.filter[["percent.mt"]] <- PercentageFeatureSet(AthRootSample.exp.filter, pattern = "^ATMG")
AthRootSample.exp.filter[["percent.chl"]] <- PercentageFeatureSet(AthRootSample.exp.filter, pattern = "^ATCG")

pdf("AthRoot.step02.01.VlnPlot.pdf", height=6, width=8)
VlnPlot(AthRootSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

AthRootSample.exp.filter <- subset(AthRootSample.exp.filter, subset = nFeature_RNA > 200 & nFeature_RNA < 9600 & percent.mt < 0.07 & percent.chl < 0.1)

pdf("AthRoot.step02.01.VlnPlot.filter.pdf", height=6, width=8)
VlnPlot(AthRootSample.exp.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"), ncol = 4)
dev.off()

AthRootSample.exp.filter <- NormalizeData(AthRootSample.exp.filter, normalization.method = "LogNormalize", scale.factor = 150000)

AthRootSample.exp.filter <- FindVariableFeatures(AthRootSample.exp.filter, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(AthRootSample.exp.filter), 10)

pdf("AthRoot.step02.02.VariableFeaturePlot.pdf", height=5, width=12)
VariableFeaturePlot(AthRootSample.exp.filter) + LabelPoints(plot = VariableFeaturePlot(AthRootSample.exp.filter), points = top10, repel = TRUE)
dev.off()

all.genes <- rownames(AthRootSample.exp.filter)
AthRootSample.exp.filter <- ScaleData(AthRootSample.exp.filter, features = all.genes)

AthRootSample.exp.filter <- RunPCA(AthRootSample.exp.filter, npcs = 50, features = VariableFeatures(object = AthRootSample.exp.filter))

pdf("AthRoot.step02.03.VizDimLoadings.pdf", height=6, width=8)
VizDimLoadings(AthRootSample.exp.filter, dims = 1:2, reduction = "pca")
dev.off()

pdf("AthRoot.step02.04.DimPlot.pdf", height=6, width=6)
DimPlot(AthRootSample.exp.filter, reduction = "pca")
dev.off()

pdf("AthRoot.step02.05.DimHeatmap.pdf", height=16, width=8)
DimHeatmap(AthRootSample.exp.filter, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

AthRootSample.exp.filter <- JackStraw(AthRootSample.exp.filter, num.replicate = 100, dims = 50)
AthRootSample.exp.filter <- ScoreJackStraw(AthRootSample.exp.filter, dims = 1:30)

pdf("AthRoot.step02.06.JackStrawPlot.pdf", height=6, width=8)
JackStrawPlot(AthRootSample.exp.filter, dims = 1:30)
dev.off()

pdf("AthRoot.step02.07.ElbowPlot.pdf", height=5, width=5)
ElbowPlot(AthRootSample.exp.filter, ndims = 30)
dev.off()

AthRootSample.exp.filter <- FindNeighbors(AthRootSample.exp.filter, dims = 1:30)

AthRootSample.exp.filter <- FindClusters(AthRootSample.exp.filter, resolution = 0.6)

AthRootSample.exp.filter.UMAP <- RunUMAP(AthRootSample.exp.filter, dims = 1:30)
AthRootSample.exp.filter.tSNE <- RunTSNE(AthRootSample.exp.filter, dims = 1:30)

pdf("AthRoot.step02.08.DimPlot.UMAP.pdf", height=3.5, width=5)
DimPlot(AthRootSample.exp.filter.UMAP, reduction = "umap")
dev.off()

pdf("AthRoot.step02.09.DimPlot.tSNE.pdf", height=3.5, width=5)
DimPlot(AthRootSample.exp.filter.tSNE, reduction = "tsne")
dev.off()

color_list<-hue_pal()(14)
new_identity_TSNE<-AthRootSample.exp.filter.tSNE
new_identity_TSNE<-RenameIdents(new_identity_TSNE, "0"="3+6", "1"="1", "2"="0", "3"="4", "4"="2", "5"="7", "6"="5", "7"="9", "8"="10", "9"="11", "10"="12", "11"="8", "12"="13")
pdf("AthRoot.step02.09.DimPlot.tSNE.select.pdf", height=4, width=4.9)
DimPlot(new_identity_TSNE, reduction = "tsne", cols=c(color_list[4],color_list[2],color_list[1],color_list[5],color_list[3],color_list[8],color_list[6],color_list[10],color_list[11],color_list[12],color_list[13],color_list[9],color_list[14]))
dev.off()

#      1 OriCluster      NewCluster	OriCluster	NewCluster
#     19 A14     B13	A14	B13
#     63 A09     B12	A09	B12
#    136 A13     B11	A13	B11
#    143 A12     B10	A12	B10
#    153 A11     B09	A11	B09
#    155 A10     B08	A10	B08
#    199 A06     B07	A06	B07
#    227 A08     B06	A08	B06
#    249 A03     B05	A03	B05
#    271 A05     B04	A05	B04
#    351 A01     B03	A01	B03
#    357 A02     B02	A02	B02
#    304 A04     B01	A04	B01
#    108 A07     B01	A07	B01

AthRootSample.exp.filter.tSNE.markers <- FindAllMarkers(AthRootSample.exp.filter.tSNE, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
AthRootSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = myAUC)

top10 <- AthRootSample.exp.filter.tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)

pdf("AthRoot.step02.10.DoHeatmap.pdf", height=8, width=8)
DoHeatmap(AthRootSample.exp.filter.tSNE, features = top10$gene) + NoLegend()
dev.off()

saveRDS(AthRootSample.exp.filter, file = "AthRoot_step02.rds")

write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="0"), "AthRoot.step02.Ath01_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="1"), "AthRoot.step02.Ath02_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="2"), "AthRoot.step02.Ath03_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="3"), "AthRoot.step02.Ath04_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="4"), "AthRoot.step02.Ath05_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="5"), "AthRoot.step02.Ath06_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="6"), "AthRoot.step02.Ath07_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="7"), "AthRoot.step02.Ath08_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="8"), "AthRoot.step02.Ath09_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="9"), "AthRoot.step02.Ath10_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="10"), "AthRoot.step02.Ath11_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="11"), "AthRoot.step02.Ath12_ID.txt", sep="\t", quote = FALSE)
write.table(WhichCells(AthRootSample.exp.filter.tSNE, idents="12"), "AthRoot.step02.Ath13_ID.txt", sep="\t", quote = FALSE)

pdf("AthRoot.step02.11.FeaturePlot.select.pdf", height=6.3, width=10)
FeaturePlot(AthRootSample.exp.filter.tSNE, features = c("AT4G36710","AT1G07640","AT1G68810","AT2G37090","AT5G57620","AT3G11550","AT2G34910","AT1G27740","AT2G37260","AT5G18840","AT5G17520","AT1G28290"))
dev.off()
