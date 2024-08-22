# Data processing from Thomas et al. study "Single-cell transcriptomic analyses reveal
# distinct immune cell contributions to epithelial barrier dysfunction in checkpoint
# inhibitor colitis" Nature Medicine, 2024
# doi: https://doi.org/10.1038/s41591-024-02895-x

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


categories_blood <- c(
    "cell", "barcode", "donor", "case", "class", "class_short",
    "ir_colitis", "drug", "type_of_cpi", "type_of_it",
    "UMAP1", "UMAP2", "leiden", "cluster", "n_counts", "n_features", "mito_counts", "mito_pct",
    "has_tcr", "TRAV", "TRAV2", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBV2", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt",
    "gender", "cancer_type", "metastatic_disease", "other_cancer_treatments"
)

metanames_blood <- c(
    "Cell_ID", "Barcodes", "Patient", "Case", "Class", "Class2",
    "irColitis", "Drug", "DrugName", "DrugType",
    "UMAP1", "UMAP2", "Leiden", "Cluster", "n_Counts", "n_Features", "mito_Counts", "mito_pct",
    "HasTCR", "TRAV", "TRAV2", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBV2", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt",
    "Gender", "CancerType", "MetastaticDisease", "OtherTreatments"
)

categories_tissue <- c(
    "cell", "barcode", "donor", "case", "class", "class_short", "drug",
    "UMAP1", "UMAP2", "leiden", "cluster", "n_counts", "n_features", "mito_counts", "mito_pct",
    "has_tcr", "TRAV", "TRAV2", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBV2", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt"
)

metanames_tissue <- c(
    "Cell_ID", "Barcodes", "Patient", "Case", "Class", "Class2", "Drug",
    "UMAP1", "UMAP2", "Leiden", "Cluster", "n_Counts", "n_Features", "mito_Counts", "mito_pct",
    "HasTCR", "TRAV", "TRAV2", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBV2", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt"
)

replace_patient <- c(
    "SIC_100" = "C1", "SIC_43" = "C10", "SIC_71" = "C11", "SIC_89" = "C12", "SIC_97" = "C13", "SIC_76" = "C14",
    "SIC_121" = "C2", "SIC_126" = "C3", "SIC_134" = "C4", "SIC_140" = "C5", "SIC_141" = "C6", "SIC_32" = "C7",
    "SIC_36" = "C8", "SIC_40" = "C9", "SIC_13" = "HC1", "SIC_186" = "HC2", "SIC_187" = "HC3", "SIC_188" = "HC4",
    "SIC_196" = "HC5", "MC_1" = "HC6", "MC_2" = "HC7", "MC_9" = "HC8", "SIC_109" = "IHC1", "SIC_172" = "IHC2",
    "SIC_19" = "IHC3", "SIC_31" = "IHC4", "SIC_53" = "IHC5", "SIC_94" = "IHC6", "SIC_33" = "IHC7", "SIC_132" = "IHC8"
)

blood_cd8_clusters <- c(
    "1" = "CX3CR1 FGFBP2 ZNF683", "2" = "HLA-DRA GZMK", "3" = "IL7R SELL", "4" = "GDT KLRC3",
    "5" = "CX3CR1 FGFBP2", "6" = "NK FCER1G", "7" = "CX3CR1 FGFBP2 TPM4", "8" = "HLA-DRA GZMK ZNF683",
    "9" = "CX3CR1 FGFBP2 NEAT1", "10" = "MAIT", "11" = "IL7R SELL CCR7", "12" = "Cycling",
    "13" = "CX3CR1 FGFBP2 TRBV12−4", "14" = "NK NCAM1/GDT/NKT"
)

tissue_cd8_clusters <- c(
    "1" = "ITGAE KIR2DL4", "2" = "ITGAE GZMB", "3" = "ITGB2 GZMK", "4" = "Cycling", "5" = "ITGAE HAVCR2",
    "6" = "ITGAE TCF7", "7" = "ITGAE GZMB", "8" = "GDT and NK", "9" = "ITGAE CD69", "10" = "CD8/4 FOXP3",
    "11" = "ITGB2 CX3CR1", "12" = "T/B cell doublet"
)

blood_cd4_clusters <- c(
    "1" = "CD4 T CD45RA CCR7", "2" = "CD4 T SELL IL7R", "3" = "CD4 T HLA-DRA PDCD1", "4" = "CD4 T LGALS3", "5" = "CD4 T GZMK CD69",
    "6" = "CD4 T", "7" = "CD4 T HLA-DRA FOXP3"
)

tissue_cd4_clusters <- c(
    "1" = "Naive/Tcm", "2" = "CT4 T ITGA1 CCL5", "3" = "CD4 T TCF7 FOS", "4" = "CT4 T ITGA1 IL23R",
    "5" = "Treg SELL", "6" = "Tfh", "7" = "Th1/Th17", "8" = "Treg KLRB1", "9" = "Treg TNFRSF18", "10" = "CD4 T Cycling"
)

# list files
files <- list.files(path = paste(getwd(), "/data/thomas/", sep = ""), pattern = "(cd4|cd8).*\\.h5ad$", ignore.case = TRUE)
# files <- "GSE206298_ircolitis-blood-b.h5ad"
meta.list <- list()
matrix.list <- list()

for (file in files) {
    # load h5ad file
    ad <- read_h5ad(paste(getwd(), "/data/thomas/", file, sep = ""))
    # extract matrix and metadata
    matrix <- t(ad$X)
    if (str_detect(file, pattern = "blood")) {
        categories <- categories_blood
        metanames <- metanames_blood
    } else if (str_detect(file, pattern = "tissue")) {
        categories <- categories_tissue
        metanames <- metanames_tissue
    }
    metadata <- ad$obs[, categories]
    names(metadata) <- metanames
    cat("\nMatrix barcodes and metadata barcodes for", file, "match:", identical(colnames(matrix), metadata$Cell_ID), "\n")
    # rename genes
    gene_names <- rownames(matrix)
    renamed_genes <- sub(".*\\|", "", gene_names)
    rownames(matrix) <- renamed_genes
    # set NAs
    metadata$Barcodes <- sub(".*\\|", "", metadata$Cell_ID)
    metadata$Patient2 <- replace_patient[match(metadata$Patient, names(replace_patient))]
    if (file == "GSE206298_ircolitis-blood-cd4.h5ad") {
        metadata$CellType <- blood_cd4_clusters[match(metadata$Leiden, names(blood_cd4_clusters))]
        metadata$Tissue <- "Blood"
        metadata$Tcell <- "CD4"
        metadata$Sample_ID <- "CD4_Blood"
    } else if (file == "GSE206298_ircolitis-blood-cd8.h5ad") {
        metadata$CellType <- blood_cd8_clusters[match(metadata$Leiden, names(blood_cd8_clusters))]
        metadata$Tissue <- "Blood"
        metadata$Tcell <- "CD8"
        metadata$Sample_ID <- "CD8_Blood"
    } else if (file == "GSE206299_ircolitis-tissue-cd4.h5ad") {
        metadata$CellType <- tissue_cd4_clusters[match(metadata$Leiden, names(tissue_cd4_clusters))]
        metadata$Tissue <- "Colon"
        metadata$Tcell <- "CD4"
        metadata$Sample_ID <- "CD4_Colon"
    } else if (file == "GSE206299_ircolitis-tissue-cd8.h5ad") {
        metadata$CellType <- tissue_cd8_clusters[match(metadata$Leiden, names(tissue_cd8_clusters))]
        metadata$Tissue <- "Colon"
        metadata$Tcell <- "CD8"
        metadata$Sample_ID <- "CD8_Colon"
    }
    # rename columns data
    metadata$Class <- gsub("Screening colonoscopy", "Healthy", metadata$Class)
    metadata$Class2 <- gsub("SC", "H", metadata$Class2)
    if (str_detect(file, pattern = "blood")) {
        metadata$irColitis <- gsub("YES", "yes", metadata$irColitis)
        metadata$irColitis <- gsub("NO", "no", metadata$irColitis)
    }
    metadata$HasTCR <- gsub("TRUE", "yes", metadata$HasTCR)
    metadata$HasTCR <- gsub("FALSE", "no", metadata$HasTCR)
    # set as factor
    metadata$Leiden <- as.factor(metadata$Leiden)
    metadata$CellType <- as.factor(metadata$CellType)
    metadata$Class <- as.factor(metadata$Class)
    metadata$Class2 <- as.factor(metadata$Class2)
    # reorder
    metadata <- metadata %>% relocate(c(CellType, Tcell), .after = Leiden)
    metadata <- metadata %>% relocate(c(Patient2, Tissue, Sample_ID), .after = Patient)
    # rename cell ids
    cell_ids <- metadata$Cell_ID
    new_ids <- paste("Thomas", metadata$Patient2, metadata$Sample_ID, metadata$Barcodes, sep = "_")
    new_ids[which(duplicated(new_ids, fromLast = TRUE))] <- paste(new_ids[which(duplicated(new_ids, fromLast = TRUE))], "_2", sep = "")
    metadata$Cell_ID <- new_ids
    rownames(metadata) <- new_ids
    colnames(matrix) <- new_ids
    # store files in list
    meta.list[[file]] <- metadata
    matrix.list[[file]] <- matrix
    # replicate plots
    ggplot(metadata, aes(x = UMAP1, y = UMAP2, color = CellType)) +
        geom_point()
    p1 <- ggplot(metadata, aes(x = n_Counts, y = n_Features)) +
        geom_point(color = "blue")
    p2 <- ggplot(metadata, aes(x = n_Counts, y = mito_pct)) +
        geom_point(color = "blue")
    p3 <- ggplot(metadata, aes(x = n_Features, y = mito_pct)) +
        geom_point(color = "blue")
    plots <- list(p1, p2, p3)
    plot_name <- gsub("^[^_]*_|\\.h5ad$", "", file)
    grid_plots <- arrangeGrob(grobs = plots, ncol = 3, nrow = 1, top = textGrob(plot_name, gp = gpar(fontsize = 16, fontface = "bold")))
    ggsave(paste(plot_name, ".pdf", sep = ""), grid_plots, width = 12, height = 6)
    rm(grid_plots, matrix, metadata, p1, p2, p3, plots, ad, cell_ids, gene_names, new_ids, renamed_genes, plot_name, metanames, categories)
}

# create unique metadata
metadata <- bind_rows(meta.list)
# add patient info to colon samples
missing_meta_cols <- setdiff(names(meta.list[[1]]), names(meta.list[[3]]))
# group info for single patients
grouped_meta <- metadata %>%
    select(Patient2, all_of(missing_meta_cols)) %>%
    group_by(Patient2) %>%
    summarise(across(everything(), ~ .[1]))
# add info for each patient and each column
for (pat in grouped_meta$Patient2) {
    for (col in missing_meta_cols) {
        metadata[[col]][which(metadata$Patient2 == pat)] <- grouped_meta[[col]][which(grouped_meta$Patient2 == pat)]
    }
}

pat_meta <- read_xlsx(path = paste(getwd(), "/data/thomas/Supplementary Tables.xlsx", sep = ""), sheet = "Patients Metadata")
for (col in colnames(pat_meta)[-1]) {
    replace_vector <- setNames(pat_meta[[col]], pat_meta$Patient2)
    metadata[[col]] <- replace_vector[match(metadata$Patient2, names(replace_vector))]
}

save(metadata, meta.list, matrix.list, file = "thomas_processed_data.RData")
write.csv(metadata, file = "test_thomas_metadata_all.csv")

# num_cells_with_tcr = length(which(x$HasTCR == 'yes'))
# num_cells_with_matched_tcr = length(which(x$PairedTCR == 'yes'))
