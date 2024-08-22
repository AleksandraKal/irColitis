setwd(getwd())

load("thomas_processed_data.RData")
identical(rownames(matrix.list[[3]]), rownames(matrix.list[[4]]))
tissue.merged <- cbind(matrix.list[[3]], matrix.list[[4]])
tissue.merged <- tissue.merged[-which(duplicated(rownames(tissue.merged))), ]
tissue.meta <- metadata[which(metadata$Tissue == "Colon"), ]
# Run Seurat pipeline
seuTissue <- CreateSeuratObject(counts = tissue.merged, meta.data = tissue.meta)
