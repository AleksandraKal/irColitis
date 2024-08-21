library(dplyr)

df1 <- read.csv("thomas_metadata_all.csv")
df2 <- read.csv("luoma_metadata.csv")
missing_columns <- setdiff(names(df1), names(df2))
for (col in missing_columns) {
    df2[[col]] <- NA
}
df2 <- df2[, names(df1)]
combined_df <- rbind(df1, df2)

write.csv(combined_df, file = "combined_file.csv", row.names = FALSE)
