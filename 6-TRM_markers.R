rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

##############
load("CD8T.rda")
table(liver$seurat_clusters2)


#AddModuleScore
DefaultAssay(liver) <- "RNA"
PL_features <- list(gene_list)

PLscore <- AddModuleScore(liver,
                          features = PL_features,
                          ctrl = 100,
                          name = "PL_Features")
colnames(PLscore@meta.data)
colnames(PLscore@meta.data)[16] <- 'TRM_Score' 

VlnPlot(PLscore,features = 'TRM_Score', 
        pt.size = 0, adjust = 2,group.by = "seurat_clusters2")

data<- FetchData(PLscore,vars = c("seurat_clusters2","TRM_Score"))
my_comparisons=list(c("0","1"))

ggplot(data, aes(x = seurat_clusters2, y = TRM_Score )) +
  geom_violin(aes(col = seurat_clusters2)) +
  geom_boxplot(aes(col = seurat_clusters2),width = 0.4)+
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")