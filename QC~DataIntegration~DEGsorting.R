library(dplyr)
library(patchwork)
library(Seurat)

#Quality control
Control <- CreateSeuratObject(counts = Control.data, project = "Control", min.cells = 3)
Gem <- CreateSeuratObject(counts = Gem.data, project = "Gem", min.cells = 3)
ATG <- CreateSeuratObject(counts = ATG.data, project = "ATG", min.cells = 3)
Control[["percent.mt"]] <- PercentageFeatureSet(Control, pattern = "^Mt-", assay = 'RNA')
VlnPlot(Control, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size=1, ncol = 3)
Control <- subset(Control, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)
Gem[["percent.mt"]] <- PercentageFeatureSet(Gem, pattern = "^Mt-", assay = 'RNA')
VlnPlot(Gem, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size=1, ncol = 3)
Gem <- subset(Gem, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)
ATG[["percent.mt"]] <- PercentageFeatureSet(ATG, pattern = "^Mt-", assay = 'RNA')
VlnPlot(ATG, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size=1, ncol = 3)
ATG <- subset(ATG, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)

#Data Merging, Normalization, Integration with barch correction
Pan02 <- merge(x=Control, y = c(Gem, ATG), add.cell.ids = c("Control","Gem","ATG"), project = "Pan02")
Pan02.list <- SplitObject(Pan02,split.by = "orig.ident")
Pan02.list <- snubh_GC.list[c("Control","Gem","ATG")]
for (i in 1:length(x = Pan02.list)) {
  Pan02.list[[i]] <- NormalizeData(object = Pan02.list[[i]], verbose = FALSE)
  Pan02.list[[i]] <- FindVariableFeatures(object = Pan02.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}
reference.list <- Pan02.list[c("Control","Gem","ATG")]
Pan02.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
Pan02_integrated <- IntegrateData(anchorset = Pan02.anchors, dims = 1:30)
DefaultAssay(Pan02_integrated) <- "integrated"

#Clustering
all.genes <- rownames(Pan02_integrated)
Pan02_integrated <- ScaleData(Pan02_integrated, features = all.genes)
Pan02_integrated <- RunPCA(Pan02_integrated, features = VariableFeatures(object =Pan02_integrated))
print(Pan02_integrated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Pan02_integrated, dims = 1:2, reduction = "pca")
DimPlot(Pan02_integrated, reduction = "pca")
DimHeatmap(Pan02_integrated, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(Pan02_integrated, ndims = 100)

Pan02_integrated <- FindNeighbors(Pan02_integrated, dims = 1:47)
Pan02_integrated <- FindClusters(Pan02_integrated, resolution =0.5)
Pan02_integrated <- RunUMAP(Pan02_integrated, dims = 1:47,umap.method='umap-learn', metric='correlation')  
DimPlot(Pan02_integrated,reduction = "umap", pt.size = 0.8, label=TRUE)

#DEG sorting
DefaultAssay(Pan02_CAF) <- "RNA"
CAF1.markers <- FindMarkers(Pan02_CAF, ident.1 = "CAF-1", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(CAF1.markers,file="CAF1.csv")
CAF2.markers <- FindMarkers(Pan02_CAF, ident.1 = "CAF-2", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(CAF2.markers,file="CAF2.csv")
CAF3.markers <- FindMarkers(Pan02_CAF, ident.1 = "CAF-3", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(CAF3.markers,file="CAF3.csv")
CAF4.markers <- FindMarkers(Pan02_CAF, ident.1 = "CAF-4", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(CAF4.markers,file="CAF4.csv")
CAF5.markers <- FindMarkers(Pan02_CAF, ident.1 = "CAF-5", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(CAF5.markers,file="CAF5.csv")
CAF6.markers <- FindMarkers(Pan02_CAF, ident.1 = "CAF-6", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(CAF6.markers,file="CAF6.csv")
  #In .csv file, 'diff.pct' was calculated by pct.1 minus pct.2
