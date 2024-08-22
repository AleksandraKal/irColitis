drop <- c("leiden0.5", "leiden0.644", "leiden0.789", "leiden0.933", "leiden1.08", "leiden1.22", "leiden1.37", "leiden1.51", "leiden1.66", "leiden1.8", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

intial_subset <- c("cell", "n_counts", "n_features", "mito_counts", "mito_pct", "donor_batch", "donor", "etiology_of_symptoms_for_controls", "patient_on_steroids", "case", "class", "confirmed_irae", "suspected_irae", "checkpoint_colitis", "checkpoint_hepatitis", "pre_immunotherapy", "gender", "cancer_type", "suspected_ir_ae_at_time_of_evaluation", "lung_cancer_subtype", "metastatic_disease", "met_location", "type_of_cpi", "type_of_it", "other_cancer_treatments", "race", "ethnicity", "ir_colitis", "class_short", "drug", "barcode", "TRAV", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt", "TRAV2", "TRBV2", "has_tcr", "UMAP1", "UMAP2", "leiden", "cluster")

subset_blood <- c(
    "cell", "barcode", "donor", "case", "class", "class_short",
    "ir_colitis", "drug", "type_of_cpi", "type_of_it",
    "UMAP1", "UMAP2", "leiden", "cluster", "n_counts", "n_features", "mito_counts", "mito_pct",
    "has_tcr", "TRAV", "TRAV2", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBV2", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt",
    "gender", "cancer_type", "metastatic_disease", "other_cancer_treatments"
)

# a and b parts: https://www.ncbi.nlm.nih.gov/books/NBK27145/

rename_blood <- c(
    "Cell_ID", "Barcodes", "Patient", "Case", "Class", "Class2",
    "irColitis", "Drug", "DrugName", "DrugType",
    "UMAP1", "UMAP2", "Leiden", "Cluster", "n_Counts", "n_Features", "mito_Counts", "mito_pct",
    "HasTCR", "TRAV", "TRAV2", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBV2", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt",
    "Gender", "CancerType", "MetastaticDisease", "OtherTreatments"
)

subset_tissue <- c(
    "cell", "barcode", "donor", "case", "class", "class_short", "drug",
    "UMAP1", "UMAP2", "leiden", "cluster", "n_counts", "n_features", "mito_counts", "mito_pct",
    "has_tcr", "TRAV", "TRAV2", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBV2", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt"
)

rename_tissue <- c(
    "Cell_ID", "Barcodes", "Patient", "Case", "Class", "Class2", "Drug",
    "UMAP1", "UMAP2", "Leiden", "Cluster", "n_Counts", "n_Features", "mito_Counts", "mito_pct",
    "HasTCR", "TRAV", "TRAV2", "TRAJ", "TRAC", "TRA_cdr3", "TRA_cdr3_nt", "TRBV", "TRBV2", "TRBD", "TRBJ", "TRBC", "TRB_cdr3", "TRB_cdr3_nt"
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

file_mapping <- list(
    "GSE206298_ircolitis-blood-cd4.h5ad" = list(
        CellType = blood_cd4_clusters,
        Tissue = "Blood",
        Tcell = "CD4",
        Sample_ID = "CD4_Blood"
    ),
    "GSE206298_ircolitis-blood-cd8.h5ad" = list(
        CellType = blood_cd8_clusters,
        Tissue = "Blood",
        Tcell = "CD8",
        Sample_ID = "CD8_Blood"
    ),
    "GSE206299_ircolitis-tissue-cd4.h5ad" = list(
        CellType = tissue_cd4_clusters,
        Tissue = "Colon",
        Tcell = "CD4",
        Sample_ID = "CD4_Colon"
    ),
    "GSE206299_ircolitis-tissue-cd8.h5ad" = list(
        CellType = tissue_cd8_clusters,
        Tissue = "Colon",
        Tcell = "CD8",
        Sample_ID = "CD8_Colon"
    )
)

replace_patient <- c(
    "SIC_100" = "C1", "SIC_43" = "C10", "SIC_71" = "C11", "SIC_89" = "C12", "SIC_97" = "C13", "SIC_76" = "C14",
    "SIC_121" = "C2", "SIC_126" = "C3", "SIC_134" = "C4", "SIC_140" = "C5", "SIC_141" = "C6", "SIC_32" = "C7",
    "SIC_36" = "C8", "SIC_40" = "C9", "SIC_13" = "HC1", "SIC_186" = "HC2", "SIC_187" = "HC3", "SIC_188" = "HC4",
    "SIC_196" = "HC5", "MC_1" = "HC6", "MC_2" = "HC7", "MC_9" = "HC8", "SIC_109" = "IHC1", "SIC_172" = "IHC2",
    "SIC_19" = "IHC3", "SIC_31" = "IHC4", "SIC_53" = "IHC5", "SIC_94" = "IHC6", "SIC_33" = "IHC7", "SIC_132" = "IHC8"
)
