# cd3 - coreceptor on T cells
# cd45 - expressed on white blood cells!
library(pacman)
p_load(Seurat, tidyverse, readxl, ggplot2, dplyr, stringr, gridExtra, grid, Matrix, anndata, purrr, MASS, dplyr, ggplot2)


## adding in the Supplemntary tables 1,2 and 3 + clusters for annotations
constPath <- paste(getwd(), "/luoma/constants.r", sep = "")
source(constPath)


filePath <- paste(getwd(), "/data/luoma/", sep = "")
folders <- list.files(filePath, pattern = "CD3")
# filter out annotations
folders <- setdiff(folders, "CD3.cell.annotation.txt")

seurat_objects <- list()
for (file in folders) {
    # print(file)
    inital <- Read10X(data.dir = paste(filePath, file, sep = ""))
    seu <- CreateSeuratObject(counts = intial, min.features = 200)
    patient <- str_extract(file, "^[A-Z]{1,2}[0-9]+")
    seu@meta.data$Patient <- patient
    patient_info <- patient_metadata[[patient]]
    for (info_name in names(patient_info)) {
        seu@meta.data[[info_name]] <- patient_info[[info_name]]
        # print(c(info_name, patient_info[[info_name]]))
    }
    count_matrix <- GetAssayData(seu, layer = "counts")
    rownames()
    seurat_objects[[file]] <- seu
}

combined_seurat <- merge(seurat_objects[[1]], y = seurat_objects[-1])
metadata <- seurat_obj@meta.data
write.csv(metadata, file = "luoma_metadata.csv", row.names = TRUE)
# path <- "/Users/aleksandrakalinic/Documents/uni/24/24T2/Thesis a/final/irColitis/data/luoma/"
# folders <- list.files(C1Path, pattern = "CD")

# intial <- Read10X(data.dir = paste(path, "C1-CD3", sep = ""))
# seu <- CreateSeuratObject(counts = intial, min.features = 200)

# ANNOTATIONS:
# annotations <- paste(getwd(), "/data/luoma/", sep = "")
# folders <- list.files(annotations, pattern = "^CD(3|4|8|45)\\.cell\\.annotation\\.txt$")
# list_of_dfs <- list()
# for (file in folders) {
#     cd_annotations <- read.table(paste(annotations, file, sep = ""), header = TRUE, sep = "")
#     cd_annotations$Barcode <- rownames(cd_annotations)
#     cd_annotations$CellType_ID <- factor(cd_annotations$CellType_ID, levels = unique(cd_annotations$CellType_ID))
#     cd_find <- cd_clusters[[file]]
#     cd_annotations$ClusterName <- cd_find[match(cd_annotations$CellType_ID, names(cd_find))]
#     list_of_dfs <- append(list_of_dfs, list(cd_annotations))
# }

# combined_df <- do.call(rbind, list_of_dfs)
# matrices <- list()
# for (file in folders) {
#     filePath <- paste(dataPath, file, sep = "")
#     print(filePath)
#     matrices[[file]] <- Read10X(data.dir = filePath)
# }

for (patient in names(patient_metadata)) {
    print(patient)
    print("=======================================")
    # add to metadata
    # seu@metadata$patient = patient
    patient_info <- patient_metadata[[patient]]
    for (info_name in names(patient_info)) {
        # seu@metadata$info_name=patient_info[[info_name]]
        print(c(info_name, patient_info[[info_name]]))
    }
}
