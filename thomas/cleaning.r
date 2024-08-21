library(pacman)

p_load(Seurat, tidyverse, readxl, ggplot2, dplyr, stringr, gridExtra, grid, Matrix, anndata)

# constPath <- paste(getwd(), "/thomas/constants.r", sep = "")
# retrive the constants
# source(constPath)

# TODO: make this more constant
path <- paste(getwd(), "/irColitis/data/thomas/", sep = "")
files <- list.files(path, pattern = "(cd4|cd8).*\\.h5ad$", ignore.case = TRUE)
files <- "GSE206299_ircolitis-tissue-cd4.h5ad"
# files <- "GSE206298_ircolitis-blood-cd4.h5ad"

meta.list <- list()
matrix.list <- list()

for (file in files) {
    # TODO: fix to a nicer format
    ad <- read_h5ad(paste(getwd(), "/data/thomas/", file, sep = ""))
    matrix <- t(ad$X)
}

# 1. change the rownames to basepair-Patient
