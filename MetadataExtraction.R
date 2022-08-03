library(dplyr)
library(Seurat)

#Extract filtered cell IDs
write.csv(Cells(Pan02_CAF), file="cellID_obs.csv", row.names=FALSE)

#Extract UMAP coordinates
write.csv(Embeddings(Pan02_CAF), reduction="umap"), file="cell_embeddings.csv")

#Extract cluster info
write.csv(Pan02_CAF@meta.data$seurat_clusters, file="clusters.csv")
