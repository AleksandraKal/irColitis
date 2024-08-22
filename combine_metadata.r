library(dplyr)

df1 <- read.csv("thomas_metadata_all.csv")
df2 <- read.csv("luoma_metadata.csv")

# df1 <- df1[order(colnames(df1)), ]
# df2 <- df2[order(df2$column_name), ]

# change thomas data set
df1$study <- "Thomas"
df1 <- subset(df1, select = -Cell_ID)
colnames(df2)[colnames(df2) == "OtherTreatments"] <- "Treatment"
colnames(df2)[colnames(df2) == "Gender"] <- "Sex"

# rename to same as thomas
colnames(df2)[colnames(df2) == "nCount_RNA"] <- "n_Counts"
colnames(df2)[colnames(df2) == "nFeature_RNA"] <- "nFeature_RNA"
colnames(df2)[colnames(df2) == "UMAP_1"] <- "UMAP1"
colnames(df2)[colnames(df2) == "UMAP_2"] <- "UMAP2"
colnames(df2)[colnames(df2) == "HaveTCR"] <- "HasTCR"
df2 <- subset(df2, select = -orig.ident)
df2$Barcodes <- sub("-(C|CT|NC)?[0-9]+$", "", df2$"X")
df2$Cell_ID <- paste("Luoma", df2$Patient, sep = "_")
df2$irColitis <- ifelse(df2$Case == "irColitis", "YES", "NO")
df2$HasTCR <- ifelse(df2$HasTCR == "yes", TRUE, FALSE)
missing_columns <- setdiff(names(df1), names(df2))


combined_df <- bind_rows(df1, df2)
write.csv(combined_df, "test_combined_file.csv", row.names = FALSE)
# for (col in missing_columns) {
#     df2[[col]] <- NA
# }
# df2 <- df2[, names(df1)]



# combined_df <- rbind(df1, df2)

# write.csv(combined_df, file = "combined_file.csv", row.names = FALSE)
