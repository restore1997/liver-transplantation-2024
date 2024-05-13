rm(list = ls())

library(data.table)
library(Seurat)
library(dplyr)

load("0-raw_total.rda")

total[["percent.mt"]] <- PercentageFeatureSet(total, pattern = "^MT-")

VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,pt.size = 0,group.by = "sample",
        cols = c("#5176A0","#B2D0E7","#E88B20","#F8C481","#589E39","#A8D175","#B39420","#EECB60",
                 "#509B96","#95C6BE","#DA5554","#F3A6A6",
                 "#af6eeb", "#DEB887", "#5F9EA0", "#E1C7A6", "#ff7d4d", "#2E8B57", "#D8BFD8", 
                 "#FA8072", "#6B8E23", "#FF6347", "#40E0D0", "#EE82EE"))
total <- subset(total, subset = nFeature_RNA > 500 & percent.mt < 25)

table(total$sample)

liver.list <- c(sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9
              ,sample10,sample11,sample12,sample13,sample14,sample15,sample16,sample17)


liver.list <- lapply(X = liver.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = liver.list)

immune.anchors <- FindIntegrationAnchors(object.list = liver.list, anchor.features = features)

immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

total <- immune.combined
table(total$seurat_clusters)

DimPlot(total, reduction = "umap",label = T,raster=FALSE,
        cols = c("#D4A5D9", "#A3D5C8", "#FAA6C0", "#E1C7A6",
                 "#B7D6E1", "#F8A5B4", "#B6DAA7", "#A6D3B1",
                 "#C3D8A2", "#D5C2E0", "#C1DAB8", "#FEE497",
                 "#A5D1C9", "#D5E1A9", "#B2C6E2", "#D7B5E7",
                 "#E2A7D4", "#D8E3A5", "#FFF2A1", "#E6B5C2",
                 "#FFB3A3", "#C6A1E3", "#C4A9D6"))
DimPlot(total, reduction = "umap",group.by = "diagnosis")
DimPlot(total, reduction = "umap",group.by = "sample")
DimPlot(total, reduction = "umap",group.by = "tissue")

library(SingleR)

load("reference.rda")
refdata <- ref.data

testdata <- GetAssayData(total, slot="data")
clusters <- total@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

total@meta.data$celltype_hpca = "NA"
for(i in 1:nrow(celltype)){
  total@meta.data[which(total@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_hpca'] <- celltype$celltype[i]}

DimPlot(total, group.by="celltype_hpca", repel=T, label=T, reduction='umap',pt.size = 1.6)

table(total$celltype_hpca)

save(total,file = "1-atlas.rda")