library(pacman)
p_load(Seurat, tidyverse, readxl, ggplot2, dplyr, stringr, gridExtra, grid, Matrix, anndata, purrr, MASS, dplyr, ggplot2)

path <- "/Users/aleksandrakalinic/Documents/uni/24/24T2/Thesis a/final/irColitis/data/luoma/C1-CD3"
mat <- Read10X(data.dir = path)
seu <- CreateSeuratObject(counts = mat, min.features = 200)
count_matrix <- GetAssayData(seu, layer = "counts")
