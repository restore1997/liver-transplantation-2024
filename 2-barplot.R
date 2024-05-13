library(ggplot2)

pred <- total@meta.data[, c("orig.ident", "MM")]

print(sort(table(pred$MM)))

unique_clusters <- unique(pred$MM)
for (cluster in unique_clusters) {
  cluster_length <- sum(pred$MM == cluster)
  pred$value[pred$MM == cluster] <- cluster_length
}

print(sort(table(pred$value)))

ggplot(data = pred, mapping = aes(x = orig.ident, y = value, fill = MM)) + 
  geom_bar(position = "fill", stat = "identity") +
  theme_bw(base_size = 20) +
  scale_fill_manual(
    values = c("#EECB60", "#DA5554", "#E1C7A6", "#5F9EA0", 
               "#E88B20", "#5176A0")
  )


pred <- total@meta.data[, c("orig.ident", "MM")]

print(sort(table(pred$orig.ident)))

unique_clusters <- unique(pred$orig.ident)
for (cluster in unique_clusters) {
  cluster_length <- sum(pred$orig.ident == cluster)
  pred$value[pred$orig.ident == cluster] <- cluster_length
}

print(sort(table(pred$value)))

ggplot(data = pred, mapping = aes(x = MM, y = value, fill = orig.ident)) + 
  geom_bar(position = "fill", stat = "identity") +
  theme_bw(base_size = 20)+
  scale_fill_manual(
    values = c("#5176A0","#B2D0E7","#E88B20","#F8C481","#589E39","#A8D175","#B39420","#EECB60",
               "#509B96","#95C6BE","#DA5554","#F3A6A6",
               "#af6eeb", "#DEB887", "#5F9EA0", "#E1C7A6", "#ff7d4d", "#2E8B57", "#D8BFD8", 
               "#FA8072", "#6B8E23", "#FF6347", "#40E0D0", "#EE82EE")
  )
