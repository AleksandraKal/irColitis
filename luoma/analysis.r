library(Seurat)
library(tidyverse)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(gridExtra)
library(grid)
library(Matrix)
library(anndata)
library(purrr)
library(DoubletFinder)
library(MASS)
library(dplyr)
library(ggplot2)



# load in the seurat_objects
load("lumoa_seu2.r")

# path <- "/Users/aleksandrakalinic/Documents/uni/24/24T2/Thesis a/final/irColitis/data/luoma/C1-CD3"
# initial <- Read10X(data.dir = path)
# seu <- CreateSeuratObject(counts = initial, min.features = 200)

path <- paste(getwd(), "/data/luoma/", sep = "")
folders <- list.files(path, pattern = "CD3")
folders <- setdiff(folders, "CD3.cell.annotation.txt")
matrices <- list()
for (file in folders) {
    matrices[[file]] <- Read10X(data.dir = paste(path, file, sep = "/"))
}

# split T cells matrices
tcells <- matrices[str_detect(names(matrices), "CD3")]
tot_tcells <- sum(sapply(tcells, ncol))
doublet_rate <- (0.009 * sapply(tcells, ncol) / 1000)


# for (mat in names(seurat_objects)) {
# seu <- seurat_objects[[mat]]
mat <- "C1-CD3"
seu <- seurat_objects[["C1-CD3"]]
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = percent.mt < 11.13)

# DoubletFinder
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)

sweep.res <- paramSweep(seu, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]
pK <- as.character(pK[1])
pK <- as.numeric(pK)
homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)
nExp_poi <- round(doublet_rate[mat] * nrow(seu@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
seu <- doubletFinder(seu, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
# }


# merge later
# combined_seurat <- merge(seurat_objects[[1]], y = seurat_objects[-1])
