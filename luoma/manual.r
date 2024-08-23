library(pacman)

## TODO: REPLACE WITH DOUBLET FINDER JUST FASTER FOR NOW TO MANUALLY

p_load(Seurat, tidyverse, readxl, ggplot2, dplyr, stringr, gridExtra, grid, Matrix, anndata)
load("lumoa_seu2.r")

thresholds_list <- list(
    "C1-CD3" = list(upperRNA = 500, upperRNA = 50000, lowerFeat = 2000, upperFear = 6000),
    "C2-CD3" = list(upperRNA = 500, upperRNA = 50000, lowerFeat = 2000, upperFear = 6000),
    "C3-CD3" = list(upper = 20000, lower = 600),
    "C4-CD3" = list(upper = 20000, lower = 600),
    "C5-CD3" = list(upper = 20000, lower = 600),
    "C6-CD3" = list(upper = 20000, lower = 600),
    "C7-CD3" = list(upper = 20000, lower = 600),
    "C8-CD3" = list(upper = 20000, lower = 600),
    "CT1-CD3" = list(upper = 20000, lower = 600),
    "CT2-CD3" = list(upper = 20000, lower = 600),
    "CT3-CD3" = list(upper = 20000, lower = 600),
    "CT4-CD3" = list(upper = 20000, lower = 600),
    "CT5-CD3" = list(upper = 20000, lower = 600),
    "CT6-CD3" = list(upper = 20000, lower = 600),
    "CT7-CD3" = list(upper = 20000, lower = 600),
    "CT8-CD3" = list(upper = 20000, lower = 600),
    "NC1-CD3" = list(upper = 20000, lower = 600),
    "NC2-CD3" = list(upper = 20000, lower = 600),
    "NC3-CD3" = list(upper = 20000, lower = 600),
    "NC4-CD3" = list(upper = 20000, lower = 600),
    "NC5-CD3" = list(upper = 20000, lower = 600),
    "NC6-CD3" = list(upper = 20000, lower = 600),
)


# each patient file eg.C1-CD3 perform manual mitcounts removal
C1_CD3 <- seurat_objects[["C1-CD3"]]
VlnPlot(object = C1_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(C1_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# Filter cells based on upper thresholds - to  3703
seurat_objects[["C1-CD3"]] <- subset(C1_CD3, subset = nCount_RNA > 500 & nCount_RNA < 50000 & nFeature_RNA > 200 & nFeature_RNA < 6000)

summary(C1_CD3$nCount_RNA)
summary(C1_CD3$nFeature_RNA)
quantile(C1_CD3$nCount_RNA, probs = c(0.01, 0.99))
quantile(C1_CD3$nFeature_RNA, probs = c(0.01, 0.99))
plot(density(C1_CD3$nCount_RNA))
plot(density(C1_CD3$nFeature_RNA))


C2_CD3 <- seurat_objects[["C2-CD3"]]
# from nCountRNA - over 50,000 nCountRNA  and nFeatureRNA over 60,00 or 70,000
VlnPlot(object = C2_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(C2_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 4091 cells
seurat_objects[["C2-CD3"]] <- subset(C2_CD3, subset = nCount_RNA > 10 & nCount_RNA < 40000 & nFeature_RNA > 200 & nFeature_RNA < 6000)


C3_CD3 <- seurat_objects[["C3-CD3"]]
# from nCountRNA - over 45,500 nCountRNA  and nFeatureRNA over 60,00
VlnPlot(object = C3_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(C3_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 4383
seurat_objects[["C3-CD3"]] <- subset(C3_CD3, subset = nCount_RNA > 0 & nCount_RNA < 45500 & nFeature_RNA > 100 & nFeature_RNA < 6000)
plot(density(C3_CD3$nCount_RNA))
plot(density(C3_CD3$nFeature_RNA))


C4_CD3 <- seurat_objects[["C4-CD3"]]
# from nCountRNA - over 40,00 nCountRNA  and nFeatureRNA over 6000
VlnPlot(object = C4_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(C4_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3587
seurat_objects[["C4-CD3"]] <- subset(C4_CD3, subset = nCount_RNA > 100 & nCount_RNA < 40000 & nFeature_RNA > 0 & nFeature_RNA < 5000)
plot(density(C4_CD3$nCount_RNA))
plot(density(C4_CD3$nFeature_RNA))



# pretty bad data for this sample
C5_CD3 <- seurat_objects[["C5-CD3"]]
# from nCountRNA - over 40,00 nCountRNA  and nFeatureRNA over 5000
VlnPlot(object = C5_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(C5_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 4097
seurat_objects[["C5-CD3"]] <- subset(C5_CD3, subset = nCount_RNA > 100 & nCount_RNA < 45500 & nFeature_RNA > 0 & nFeature_RNA < 5500)
plot(density(C5_CD3$nCount_RNA))
plot(density(C5_CD3$nFeature_RNA))


C6_CD3 <- seurat_objects[["C6-CD3"]]
# from nCountRNA - over 50,00 nCountRNA  and nFeatureRNA over 5000
VlnPlot(object = C6_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(C6_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3206
seurat_objects[["C6-CD3"]] <- subset(C6_CD3, subset = nCount_RNA > 100 & nCount_RNA < 50000 & nFeature_RNA > 0 & nFeature_RNA < 5000)
plot(density(C6_CD3$nCount_RNA))
plot(density(C6_CD3$nFeature_RNA))

C7_CD3 <- seurat_objects[["C7-CD3"]]
# from nCountRNA - over 45,000 nCountRNA  and nFeatureRNA over 6000
VlnPlot(object = C7_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(C7_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3654
seurat_objects[["C7-CD3"]] <- subset(C7_CD3, subset = nCount_RNA > 100 & nCount_RNA < 45000 & nFeature_RNA > 0 & nFeature_RNA < 6000)
plot(density(C6_CD3$nCount_RNA))
plot(density(C6_CD3$nFeature_RNA))

C8_CD3 <- seurat_objects[["C8-CD3"]]
# from nCountRNA - over 70,000 nCountRNA  and nFeatureRNA over 6000
VlnPlot(object = C8_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(C8_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3705
seurat_objects[["C8-CD3"]] <- subset(C8_CD3, subset = nCount_RNA > 100 & nCount_RNA < 70000 & nFeature_RNA > 0 & nFeature_RNA < 6000)
plot(density(C6_CD3$nCount_RNA))
plot(density(C6_CD3$nFeature_RNA))

CT1_CD3 <- seurat_objects[["CT1-CD3"]]
# from nCountRNA - over 20,000 nCountRNA  and nFeatureRNA over 5000
VlnPlot(object = CT1_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(CT1_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3900
seurat_objects[["CT1_CD3"]] <- subset(CT1_CD3, subset = nCount_RNA > 100 & nCount_RNA < 20500 & nFeature_RNA > 0 & nFeature_RNA < 5000)
plot(density(CT1_CD3$nCount_RNA))
plot(density(CT1_CD3$nFeature_RNA))

CT2_CD3 <- seurat_objects[["CT2-CD3"]]
# from nCountRNA - over 20,000 nCountRNA  and nFeatureRNA over 4000
VlnPlot(object = CT2_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(CT2_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 2859
seurat_objects[["CT2_CD3"]] <- subset(CT2_CD3, subset = nCount_RNA > 100 & nCount_RNA < 20000 & nFeature_RNA > 0 & nFeature_RNA < 4000)
plot(density(CT2_CD3$nCount_RNA))
plot(density(CT2_CD3$nFeature_RNA))


# Load the CT3_CD3 Seurat object from the list
CT3_CD3 <- seurat_objects[["CT3-CD3"]]

# nCountrna: 100 - 25,000
# nFetaureRNA: 50 - 4050
# Visualize nCount_RNA and nFeature_RNA using Violin Plot and Feature Scatter Plot
VlnPlot(object = CT3_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(CT3_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")

# Filter the Seurat object based on nCount_RNA and nFeature_RNA thresholds
# 3837
seurat_objects[["CT3_CD3"]] <- subset(CT3_CD3, subset = nCount_RNA > 100 & nCount_RNA < 25000 & nFeature_RNA > 50 & nFeature_RNA < 4050)

# Plot density plots for nCount_RNA and nFeature_RNA
plot(density(CT3_CD3$nCount_RNA))
plot(density(CT3_CD3$nFeature_RNA))


# nothing filters
CT4_CD3 <- seurat_objects[["CT4-CD3"]]

CT5_CD3 <- seurat_objects[["CT5-CD3"]]


# ncountRNA: 0 - 20,000
# nfeatureRNA:0 - 4000
CT6_CD3 <- seurat_objects[["CT6-CD3"]]

CT7_CD3 <- seurat_objects[["CT7-CD3"]]
# from nCountRNA - over 20,000 nCountRNA  and nFeatureRNA over 4000
VlnPlot(object = CT7_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(CT7_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3563
seurat_objects[["CT7_CD3"]] <- subset(CT7_CD3, subset = nCount_RNA > 0 & nCount_RNA < 20000 & nFeature_RNA > 0 & nFeature_RNA < 4000)
plot(density(CT7_CD3$nCount_RNA))
plot(density(CT7_CD3$nFeature_RNA))

# 199 - 25,000
# nfeatureRNA: 0 4000
CT8_CD3 <- seurat_objects[["CT8-CD3"]]
# from nCountRNA - over 20,000 nCountRNA  and nFeatureRNA over 4000
VlnPlot(object = CT8_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(CT8_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3370
seurat_objects[["CT8_CD3"]] <- subset(CT8_CD3, subset = nCount_RNA > 200 & nCount_RNA < 25000 & nFeature_RNA > 0 & nFeature_RNA < 4000)
plot(density(CT8_CD3$nCount_RNA))
plot(density(CT8_CD3$nFeature_RNA))


NC1_CD3 <- seurat_objects[["NC1-CD3"]]
# from nCountRNA - over 30,000 nCountRNA  and nFeatureRNA over 5500
VlnPlot(object = NC1_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(NC1_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3081
seurat_objects[["NC1_CD3"]] <- subset(NC1_CD3, subset = nCount_RNA > 100 & nCount_RNA < 30000 & nFeature_RNA > 100 & nFeature_RNA < 5000)
plot(density(NC1_CD3$nCount_RNA))
plot(density(NC1_CD3$nFeature_RNA))


NC2_CD3 <- seurat_objects[["NC2-CD3"]]

NC3_CD3 <- seurat_objects[["NC3-CD3"]]
# from nCountRNA - over 20,000 nCountRNA  and nFeatureRNA 1000 over 5000
VlnPlot(object = NC3_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(NC3_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 2529
seurat_objects[["NC3_CD3"]] <- subset(NC3_CD3, subset = nCount_RNA > 0 & nCount_RNA < 20000 & nFeature_RNA > 1000 & nFeature_RNA < 5000)
plot(density(NC3_CD3$nCount_RNA))
plot(density(NC3_CD3$nFeature_RNA))


NC3_CD3 <- seurat_objects[["NC3-CD3"]]
# from nCountRNA - over 20,000 nCountRNA  and nFeatureRNA 1000 over 5000
VlnPlot(object = NC3_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(NC3_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 2529
seurat_objects[["NC3_CD3"]] <- subset(NC3_CD3, subset = nCount_RNA > 0 & nCount_RNA < 20000 & nFeature_RNA > 1000 & nFeature_RNA < 5000)
plot(density(NC3_CD3$nCount_RNA))
plot(density(NC3_CD3$nFeature_RNA))



NC4_CD3 <- seurat_objects[["NC4-CD3"]]


NC5_CD3 <- seurat_objects[["NC5-CD3"]]
# from nCountRNA - over 50,000 nCountRNA  and nFeatureRNA over 50 -4500
VlnPlot(object = NC5_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(NC5_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3563
seurat_objects[["NC5_CD3"]] <- subset(NC5_CD3, subset = nCount_RNA > 50 & nCount_RNA < 50000 & nFeature_RNA > 1000 & nFeature_RNA < 4500)
plot(density(NC5_CD3$nCount_RNA))
plot(density(NC5_CD3$nFeature_RNA))


NC6_CD3 <- seurat_objects[["NC6-CD3"]]
# from nCountRNA - over 20,000 nCountRNA  and nFeatureRNA over 100 - 4000
VlnPlot(object = NC6_CD3, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(NC6_CD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
# 3563
seurat_objects[["NC6_CD3"]] <- subset(NC6_CD3, subset = nCount_RNA > 100 & nCount_RNA < 20500 & nFeature_RNA > 1000 & nFeature_RNA < 4000)
plot(density(NC6_CD3$nCount_RNA))
plot(density(NC6_CD3$nFeature_RNA))
save(seurat_objects, file = "lumoa_seu2_filtered.r")
