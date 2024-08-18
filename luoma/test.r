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
library(MASS)
library(dplyr)
library(ggplot2)


cd3_annotations <- read.table(paste(getwd(), "/data/luoma/CD3.cell.annotation.txt", sep = ""), header = TRUE, sep = "")
cd3_annotations$Barcode <- rownames(cd3_annotations)
cd3_annotations$CellType_ID <- factor(cd3_annotations$CellType_ID, levels = unique(cd3_annotations$CellType_ID))
cd3_clusters <- c("1" = "MAIT", "2" = "Type 3 cytokines Trm", "3" = "Type 1 cytokines Trm", "4" = "Naive/CM", "5" = "5", "6" = "T follicular helper", "7" = "Treg", "8" = "Th1 effector", "9" = "Cycling", "10" = "Cytotoxic Effector", "11" = "LP Trm", "12" = "CD8 Trm", "13" = "IEL: CD8/gdT")
cd3_annotations$ClusterName <- cd3_clusters[match(cd3_annotations$CellType_ID, names(cd3_clusters))]

cd4_annotations <- read.table(paste(getwd(), "/data/luoma/CD4.cell.annotation.txt", sep = ""), header = TRUE, sep = "")
cd4_annotations$Barcode <- rownames(cd4_annotations)
cd4_annotations$CellType_ID <- factor(cd4_annotations$CellType_ID, levels = unique(cd4_annotations$CellType_ID))
cd4_clusters <- c("1" = "Th17 Trm", "2" = "Th1 Trm-A", "3" = "Th1 Trm-R", "4" = "Naive", "5" = "Central Memory", "6" = "Tfh", "7" = "Treg", "8" = "Cytotoxic", "9" = "Th1 Effector", "10" = "Cycling")
cd4_annotations$ClusterName <- cd4_clusters[match(cd4_annotations$CellType_ID, names(cd4_clusters))]

cd8_annotations <- read.table(paste(getwd(), "/data/luoma/CD8.cell.annotation.txt", sep = ""), header = TRUE, sep = "")
cd8_annotations$Barcode <- rownames(cd8_annotations)
cd8_annotations$CellType_ID <- factor(cd8_annotations$CellType_ID, levels = unique(cd8_annotations$CellType_ID))
cd8_clusters <- c("1" = "Trm - IEL", "2" = "Trm - LP1", "3" = "Trm - LP2", "4" = "MAIT", "5" = "Term effector", "6" = "CM/Naive", "7" = "Cytotoxic Effector", "8" = "Cycling")
cd8_annotations$ClusterName <- cd8_clusters[match(cd8_annotations$CellType_ID, names(cd8_clusters))]

metadata <- merge(cd4_annotations, cd8_annotations, by = "Barcode", all.x = TRUE, all.y = TRUE)
# removing the columns
metadata[, which(str_detect(names(metadata), "Group|Patient"))] <- NULL

names(metadata) <- gsub("[.]x", "_CD4ann", names(metadata))
names(metadata) <- gsub("[.]y", "_CD8ann", names(metadata))
metadata <- merge(cd3_annotations, metadata, by = "Barcode", all.x = TRUE, all.y = TRUE)
cd3_annotations$TcellType <- "CD3"
# annotate the cd3 cols if cd4 or cd8
cd3_annotations$TcellType[which(cd3_annotations$Barcode %in% cd4_annotations$Barcode)] <- "CD4"
cd3_annotations$TcellType[which(cd3_annotations$Barcode %in% cd8_annotations$Barcode)] <- "CD8"

tcr_sequences <- read.csv(file = paste(getwd(), "/data/luoma/GSE144469_TCR_filtered_contig_annotations_all.csv", sep = "")) %>%
    filter(is_cell == "True" & high_confidence == "True" & full_length == "True" & productive == "True") %>%
    dplyr::select(barcode, chain, v_gene, d_gene, j_gene, cdr3, reads)
unique(tcr_sequences$chain)

alpha <- tcr_sequences %>%
    filter(chain == "TRA") %>%
    group_by(barcode) %>%
    top_n(n = 4, wt = reads) %>%
    ungroup() %>%
    dplyr::select(-chain, -d_gene, -reads) %>%
    group_by(barcode) %>%
    mutate(row_number = row_number()) %>%
    pivot_wider(names_from = row_number, values_from = c(v_gene, j_gene, cdr3), names_glue = "{.value}.{row_number}") %>%
    ungroup() %>%
    relocate(c(v_gene.1, j_gene.1, cdr3.1, v_gene.2, j_gene.2, cdr3.2), .after = barcode)
new_names <- c("v_gene" = "TRAV", "j_gene" = "TRAJ", "cdr3" = "CDR3alpha")
names(alpha) <- str_replace_all(names(alpha), new_names)


beta <- tcr_sequences %>%
    filter(chain == "TRB") %>%
    group_by(barcode) %>%
    top_n(n = 4, wt = reads) %>%
    ungroup() %>%
    dplyr::select(-chain, -d_gene, -reads) %>%
    group_by(barcode) %>%
    mutate(row_number = row_number()) %>%
    pivot_wider(names_from = row_number, values_from = c(v_gene, j_gene, cdr3), names_glue = "{.value}.{row_number}") %>%
    ungroup() %>%
    relocate(c(v_gene.1, j_gene.1, cdr3.1, v_gene.2, j_gene.2, cdr3.2), .after = barcode)
new_names <- c("v_gene" = "TRBV", "d_gene" = "TRBD", "j_gene" = "TRBJ", "cdr3" = "CDR3beta")
names(beta) <- str_replace_all(names(beta), new_names)
