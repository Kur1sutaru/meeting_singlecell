setwd("D:/baylor/TRIGEMINAL/snijunsn")
library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(dplyr)
library(patchwork)
library(ggplot2)

# load data

snijunsn <- Read10X_h5(filename = "D:/baylor/TRIGEMINAL/snijunsn/b6tgnu.h5")
snijunsn <-LoadH5Seurat("D:/baylor/TRIGEMINAL/snijunsn/Neuron_paper_cell_label/objintegrated.h5seurat")

snijunsn <- CreateSeuratObject(counts = snijunsn, project = "snijunsn", min.cells = 3, min.features = 200)
snijunsn

dense.size <- object.size(snijunsn)
dense.size

sparse.size <- object.size(snijunsn)
sparse.size

snijunsn[["percent.mt"]] <- PercentageFeatureSet(snijunsn, pattern = "^MT-|^Mt-")

# Visualize QC metrics as a violin plot
jpeg(file="snijunsnPercentageFeatureSet.jpeg")
VlnPlot(snijunsn, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file="snijunsnFeatureScatter.jpeg")
plot1 <- FeatureScatter(snijunsn, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(snijunsn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


snijunsn <- NormalizeData(snijunsn, normalization.method = "LogNormalize", scale.factor = 10000)

snijunsn <- FindVariableFeatures(snijunsn, selection.method = "vst", nfeatures = 32000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(snijunsn), 10)
write.csv(top10,"snijunsntophighvariablegenes.csv")

top5 <- head(VariableFeatures(snijunsn), 5)
write.csv(top5,"snijunsntophighvariablegenes5.csv")

# plot variable features with and without labels
jpeg(file="snijunsnvariable_plot1.jpeg")
plot1 <- VariableFeaturePlot(snijunsn)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# plot variable features with and without labels
jpeg(file="snijunsnvariable5_plot1.jpeg")
plot1 <- VariableFeaturePlot(snijunsn)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(snijunsn)
snijunsn <- ScaleData(snijunsn, features = all.genes)

snijunsn <- RunPCA(snijunsn, features = VariableFeatures(object = snijunsn))

jpeg(file="snijunsnvizdim3load.jpeg")
VizDimLoadings(snijunsn, dims = 1:2, reduction = "pca")
dev.off()

snijunsn <- FindNeighbors(snijunsn, dims = 1:30)
snijunsn <- FindClusters(snijunsn, resolution = 0.2)
snijunsn <- RunUMAP(snijunsn, dims = 1:30)

jpeg(file="nuclei_combined0.230dims.jpeg")
DimPlot(snijunsn, group.by = "name",reduction = "umap")
dev.off()

jpeg(file="umap_lens_combined0.110dims.jpeg")
DimPlot(onhannotated,reduction = "umap")
dev.off()

snijunsn <- FindClusters(snijunsn, resolution = 0.1)
snijunsn <- RunUMAP(snijunsn, dims = 1:10)
jpeg(file="snijunsn0.310dims.jpeg")
DimPlot(snijunsn, reduction = "umap")
dev.off()

snijunsn <- FindClusters(snijunsn, resolution = 0.4)
snijunsn <- RunUMAP(snijunsn, dims = 1:10)
jpeg(file="snijunsn0.410dims.jpeg")
DimPlot(snijunsn, reduction = "umap", label=T)
dev.off()

# find all markers of cluster 2
cluster2.markers <- FindMarkers(snijunsn, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
write.csv(cluster2.markers, "snijunsnmarkers.csv")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
onh.markers <- FindAllMarkers(onhannotated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
onh.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(onh.markers, "onh_combinedclustermarkers.csv")



#Aditional steps
# Load the Seurat package and your Seurat object
library(Seurat)
my_seurat <- Read10X(data.dir = "path/to/data")
my_seurat <- CreateSeuratObject(counts = my_seurat$`Gene Expression`)

# Identify cell types using clustering and/or manual annotation
# For example:
my_seurat <- FindNeighbors(my_seurat, dims = 1:10)
my_seurat <- FindClusters(my_seurat, resolution = 0.5)
my_seurat$cell_type <- Idents(my_seurat)

# Create the cell type proportion barplot
VlnPlot(my_seurat, features = "cell_type", pt.size = 0, group.by = "cell_type")
VlnPlot(my_seurat, features = "cell_type", pt.size = 0, group.by = "cell_type", geom = "bar")

#Feature expression heatmap

jpeg(file="objB6TGNuheatmap.jpeg")
DoHeatmap(
  snijunsn,
  features = NULL,
  cells = NULL,
  group.by = "name",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
dev.off()

#Heatmap with a given gene list

counts <- GetAssyData(snijunsn, assay="RNA", slot="data")
genes <- c("Iba1","Gfap","Iba1","Atf3","cFos","Fos","Il1b","Il6","Nos2","Nox","Ccl2","Cd68","Itgam","Tac1","Calca","Trpm8","Trpv1","S100b","Gfra2",
           "Pou4F2","Gal","Cd55","Scn11a","Fxyd7","Ngfr","Nefh","Hapln4","Cbln2","Kcnab1","Sst","Il31ra","Apoe","Fabp7","Mpz","Gldn","Scn7a","Dcn","Pdgfra","Mgp","Alpl",
           "Cd74","Igfbp7","Tinagl1","Htr1f","Klf6","Klf9","Mt1","Egr1","Egr2","Cyr61","Nr4a1","Ctgf","Jag1","Lrp1")
counts <- as.matrix(counts[rownames(counts) %in% genes, ])

jpeg(file="genelistB6TGNuheatmap.jpeg")
DoHeatmap(
  counts,
  features = NULL,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)


dev.off()



saveRDS(snijunsn, file = "snijunsn.rds")


## convert seurat to h5ad
library(Seurat)
library(SeuratData)
library(SeuratDisk)

UpdateSeuratObject(snijunsn)
SaveH5Seurat(snijunsn, filename = "snijunsn.h5Seurat")
Convert("snijunsn.h5Seurat", dest = "h5ad")
Convert(snijunsn, dest = "h5ad")
