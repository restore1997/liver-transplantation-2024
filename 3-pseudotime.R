rm(list = ls())

library(monocle)
library(Seurat)
library(ggplot2)
library(dplyr)

load("CD8T.rda")
load("markers.rda")

counts.data <- as(as.matrix(liver@assays$RNA@data), 'sparseMatrix')
pheno.data <- liver@meta.data
feature.data <- data.frame(gene_short_name = row.names(counts.data), row.names = row.names(counts.data))

pd <- new("AnnotatedDataFrame", data = pheno.data)
fd <- new("AnnotatedDataFrame", data = feature.data)

HSMM <- newCellDataSet(counts.data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
print(head(pData(HSMM)))
table(HSMM$seurat_clusters)

HSMM <- reduceDimension(HSMM,
                        norm_method = 'log',
                        reduction_method = 'tSNE',
                        verbose = T)
HSMM <- clusterCells(HSMM)
#plot_clusters(HSMM)
ordering_genes <- row.names (subset(liver.markers, p_val_adj < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
save(HSMM,file = "3-pseudotime.rda")
load("3-pseudotime.rda")
HSMM <- orderCells(HSMM)

save(HSMM,file = "3-pseudotime.rda")

plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")#+  scale_color_gradient2 (low="#fcfbfd", mid="#9e9ac8", high="#3f007d")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "MM2")

plot_cell_trajectory(HSMM, color_by = 'seurat_clusters') + facet_wrap(~seurat_clusters)
plot_cell_trajectory(HSMM, color_by = 'MM2') + facet_wrap(~MM2)

HSMM$TRM <- "Non_TRM"
HSMM$TRM[HSMM$MM2%in%"TRM"] <- "TRM"
plot_cell_trajectory(HSMM, color_by = "TRM")

plot_cell_trajectory(HSMM, color_by = "Memory")
plot_cell_trajectory(HSMM, color_by = "tissue")
plot_cell_trajectory(HSMM, color_by = 'tissue') + facet_wrap(~tissue)
