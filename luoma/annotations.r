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

# for (patient in names(patient_metadata)) {
#     print(patient)
#     print("=======================================")
#     # add to metadata
#     # seu@metadata$patient = patient
#     patient_info <- patient_metadata[[patient]]
#     for (info_name in names(patient_info)) {
#         # seu@metadata$info_name=patient_info[[info_name]]
#         print(c(info_name, patient_info[[info_name]]))
#     }
# }
