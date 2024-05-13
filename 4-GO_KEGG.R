rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)

load("markers.rda")

genes <- rownames(total.markers[total.markers$cluster%in%"n",])
ego_ALL <- enrichGO(gene          =genes,
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
barplot(ego_ALL,showCategory = 10)


#KEGG
genelist <- bitr(genes, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
library(dplyr)
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
barplot(ekegg, showCategory=20)

save(ego_ALL,ekegg,file = "GO_KEGG.rda")
