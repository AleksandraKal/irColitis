# cd3 - coreceptor on T cells
# cd45 - expressed on white blood cells!
library(pacman)
p_load(Seurat, tidyverse, readxl, ggplot2, dplyr, stringr, gridExtra, grid, Matrix, anndata, purrr, MASS, dplyr, ggplot2)


## adding in the Supplemntary tables 1,2 and 3 + clusters for annotations
constPath <- paste(getwd(), "/luoma/constants.r", sep = "")
source(constPath)
# add in annotations
annPath <- paste(getwd(), "/lumoa_annotations.r", sep = "")
metadata_annotations <- readRDS(annPath)
metadata_annotations <- as.data.frame(metadata_annotations)


filePath <- paste(getwd(), "/data/luoma/", sep = "")
folders <- list.files(filePath, pattern = "CD3")
# filter out "CD3.cell" annotations
folders <- setdiff(folders, "CD3.cell.annotation.txt")
# folders <- "NC4-CD3"
# folders <- "C1-CD3"
seurat_objects <- list()

for (file in folders) {
    # print(file)
    initial <- Read10X(data.dir = paste(filePath, file, sep = ""))
    seu <- CreateSeuratObject(counts = initial, min.features = 200)
    patient <- str_extract(file, "^[A-Z]{1,2}[0-9]+")
    patient_info <- patient_metadata[[patient]]
    for (info_name in names(patient_info)) {
        seu@meta.data[[info_name]] <- patient_info[[info_name]]
        # print(c(info_name, patient_info[[info_name]]))
        if (info_name == "Class" && !is.na(patient_info[[info_name]])) {
            seu@meta.data$IrColitis <- patient_info[[info_name]] == "irColitis"
        }
    }
    count_matrix <- GetAssayData(seu, layer = "counts")

    # make similar to thomas data set
    seu@meta.data$Patient <- patient
    # rename cell names
    seu <- RenameCells(seu, new.names = gsub("1", patient, Cells(seu)))
    seu@meta.data$orig.ident <- paste("Luoma", patient, sep = "_")
    seu@meta.data$Tissue <- "Colon"

    # add in annotations, filter by patient
    # cell_names <- rownames(seu@meta.data)
    # # filtered_annotations <- metadata_annotations[metadata_annotations$Patient_ID == patient, ]
    # filtered_annotations <- metadata_annotations[metadata_annotations$Barcode %in% cell_names, ]
    # common_cells <- intersect(rownames(seu@meta.data), filtered_annotations$Barcode)
    # # subset only cells with only CD3, CD4, CD8 annotations
    # seu <- subset(seu, cells = common_cells)
    # # combine the metadata and filtered annotations togther,
    # combined_metadata <- cbind(seu@meta.data, filtered_annotations[, !names(filtered_annotations) %in% c("Barcode")])
    # seu@meta.data <- combined_metadata

    # ######################################

    # TODO: ADD IN NA in the coloumns!
    cell_names <- rownames(seu@meta.data)
    filtered_annotations <- metadata_annotations[metadata_annotations$Barcode %in% cell_names, ]
    missing_cells <- setdiff(rownames(seu@meta.data), filtered_annotations$Barcode)
    for (cell in missing_cells) {
        missing_annotationn <- missing_annotations <- data.frame(
            Barcode = cell,
            Patient_ID = patient,
            Group = mapClass2[[patient_info[["Class2"]]]],
            CellType_ID = NA,
            UMAP_1 = NA,
            UMAP_2 = NA,
            ClusterName = NA,
            CellType_ID_CD4ann = NA,
            UMAP_1_CD4ann = NA,
            UMAP_2_CD4ann = NA,
            ClusterName_CD4ann = NA,
            CellType_ID_CD8ann = NA,
            UMAP_1_CD8ann = NA,
            UMAP_2_CD8ann = NA,
            ClusterName_CD8ann = NA,
            HaveTCR = "no",
            PairedTCR = NA,
            TRAV.1 = NA,
            TRAJ.1 = NA,
            CDR3alpha.1 = NA,
            TRAV.2 = NA,
            TRAJ.2 = NA,
            CDR3alpha.2 = NA,
            TRBV.1 = NA,
            TRBJ.1 = NA,
            CDR3beta.1 = NA,
            TRBV.2 = NA,
            TRBJ.2 = NA,
            CDR3beta.2 = NA,
            TRGV.1 = NA,
            TRGJ.1 = NA,
            CDR3gamma.1 = NA,
            TRGV.2 = NA,
            TRGJ.2 = NA,
            CDR3gamma.2 = NA,
            TRGV.3 = NA,
            TRGJ.3 = NA,
            CDR3gamma.3 = NA,
            TRDV.1 = NA,
            TRDD.1 = NA,
            TRDJ.1 = NA,
            CDR3delta.1 = NA,
            TRDV.2 = NA,
            TRDD.2 = NA,
            TRDJ.2 = NA,
            CDR3delta.2 = NA,
            TRDV.3 = NA,
            TRDD.3 = NA,
            TRDJ.3 = NA,
            CDR3delta.3 = NA
        )
        filtered_annotations <- rbind(filtered_annotations, missing_annotations)
    }
    combined_metadata <- cbind(seu@meta.data, filtered_annotations[, !names(filtered_annotations) %in% c("Barcode", "Patient_ID")])
    seu@meta.data <- combined_metadata
    #########################################

    # cell_names <- colnames(seu)
    # existing_metadata <- seu@meta.data
    # filtered_annotations <- metadata_annotations[metadata_annotations$Barcode %in% cell_names, ]
    # filtered_annotations <- data.frame(
    #     matrix(NA, nrow = nrow(existing_metadata), ncol = ncol(filtered_annotations)), # NA for all initially
    #     row.names = rownames(existing_metadata) # Use cell names from Seurat object
    # )
    # combined_metadata <- cbind(seu@meta.data, filtered_annotations)
    # ######################################

    seurat_objects[[file]] <- seu
}

combined_seurat <- merge(seurat_objects[[1]], y = seurat_objects[-1])
combined_seurat@meta.data$study <- "Luoma"
metadata <- combined_seurat@meta.data
write.csv(metadata, file = "luoma_metadata.csv", row.names = TRUE)
