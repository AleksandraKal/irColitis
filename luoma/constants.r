# Metadata extracted from: https://www.cell.com/cms/10.1016/j.cell.2020.06.001/attachment/54042bfb-a29e-4845-b060-3cceccf68218/mmc1


# list = (
#     String::Treatment = {"Ipilimumab ", "Ipilimumab + Nivolumab" }
#     Num::Age =
#     String::Sex = M | F
#     String::Malignancy = "Melanoma"
#     Num::Stage =
#     String::Prior_Tx = {"Pembro", "Nivo", "None", "Pembro, Ipi", "Pembro, Durva", "Nivo, Ipi+Nivo"}
#     Time_from_last_Tx =
#     Class {"Colitis", "No Colitis", "Control"}
#     Diagnosis  = {"colitis", "Entero-colitis","enteritis" "Normal"}
#     CTCAE = 0-5
#     MES (Mayo Endoscopic Score ) = 0-5
#     RECIST_response = {PD= progressive disease, PR= partial response, SD= stable disease}
#     Overall_survival (months )= number (Overall Survival (months) from CPI initiation to)
#     Status = {"deceased", "alive"}
#     Other_irAEs = {"None"(Ici therapy and doesnt have - string), NA(healthy), "irEntero-colitis", "irEnteritis"}
#     Drug = {Anti-CTLA-4,"Anti-PD-1, Anti-CTLA-4",  "None"},
#     DrugName = {Ipilimumab (Anti-CTLA4), "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)", NA},
#     DrugType = {"Monotherapy", "Combination", NA}
# )



# list(
#     Treatment =
#     Age =
#     Sex =
#     Malignancy =
#     Stage =
#     Prior_Tx =
#     Time_from_last_Tx =
#     case <-
#     Diagnosis <-
#     CTCAE <-
#     MES <-
#     RECIST_response <-
#     Overall_survival <-
#     Status <-
# )

mapClass_short <- list(
    "IRC" = "Colitis",
    "ICI" = "No-Colitis",
    "H" = "Control"
)

patient_metadata <- list(
    "C1" = list(
        Treatment = "Ipilimumab",
        Age = 62.0,
        Sex = "F",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Pembro",
        Time_from_last_Tx = 0.5,
        Class = "irColitis",
        Diagnosis = "colitis",
        CTCAE = 2,
        MES = 1,
        RECIST_response = "PD",
        Overall_survival = 4.2,
        Status = "deceased",
        # created from assumptions - TODO: clean out to decrease depedencies
        Case = "Case",
        Class_short = "IRC",
        Other_irAEs = "None",
        Drug = "Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4)",
        DrugType = "Monotherapy"
    ),
    "C2" = list(
        Treatment = "Ipilimumab",
        Age = 71.2,
        Sex = "F",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Pembro",
        Time_from_last_Tx = 7,
        Class = "irColitis",
        Diagnosis = "colitis",
        CTCAE = 2,
        MES = 1,
        RECIST_response = "PR",
        Overall_survival = 13.0,
        Status = "deceased",
        Case = "Case",
        Class_short = "IRC",
        Other_irAEs = "None",
        Drug = "Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4)",
        DrugType = "Monotherapy"
    ),
    "C3" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 56.2,
        Sex = "M",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Pembro",
        Time_from_last_Tx = 0.5,
        Class = "irColitis",
        Diagnosis = "colitis",
        CTCAE = 3,
        MES = 2,
        RECIST_response = "PR",
        Overall_survival = 17.5,
        Status = "alive",
        Case = "Case",
        Class_short = "IRC",
        Other_irAEs = "None",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "C4" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 61.5,
        Sex = "M",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Nivo",
        Time_from_last_Tx = 0.5,
        Class = "irColitis",
        Diagnosis = "entero-colitis",
        CTCAE = 2,
        MES = 2,
        RECIST_response = "PD",
        Overall_survival = 6.1,
        Status = "deceased",
        Case = "Case",
        Class_short = "IRC",
        Other_irAEs = "irEntero-colitis",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "C5" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 72.8,
        Sex = "M",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "none",
        Time_from_last_Tx = NA,
        Class = "irColitis",
        Diagnosis = "entero-colitis",
        CTCAE = 1,
        MES = 2,
        RECIST_response = "SD",
        Overall_survival = 13.9,
        Status = "alive",
        Case = "Case",
        Class_short = "IRC",
        Other_irAEs = "irEntero-colitis",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "C6" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 78.4,
        Sex = "M",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "none",
        Time_from_last_Tx = NA,
        Class = "irColitis",
        Diagnosis = "entero-colitis",
        CTCAE = 3,
        MES = 1,
        RECIST_response = "PR",
        Overall_survival = 13.2,
        Status = "alive",
        Case = "Case",
        Class_short = "IRC",
        Other_irAEs = "irEntero-colitis",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "C7" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 61.6,
        Sex = "F",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Pembro",
        Time_from_last_Tx = 24,
        Class = "irColitis",
        Diagnosis = "colitis",
        CTCAE = 4,
        MES = 2,
        RECIST_response = "PD",
        Overall_survival = 5.7,
        Status = "deceased",
        Case = "Case",
        Class = "irColitis",
        Class_short = "IRC",
        Other_irAEs = "None",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "C8" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 49.5,
        Sex = "F",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Pembro",
        Time_from_last_Tx = 1,
        Class = "irColitis",
        Diagnosis = "colitis",
        CTCAE = 3,
        MES = 2,
        RECIST_response = "PD",
        Overall_survival = 11.1,
        Status = "deceased",
        Case = "Case",
        Class_short = "IRC",
        Other_irAEs = "None",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "NC1" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 43.6,
        Sex = "M",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "none",
        Time_from_last_Tx = NA,
        Class = "On ICI Therapy",
        Diagnosis = "enteritis",
        CTCAE = 1,
        MES = 0,
        RECIST_response = "PR",
        Overall_survival = 17.5,
        Status = "alive",
        Case = "Control",
        Class_short = "ICI",
        Other_irAEs = "irEnteritis",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "NC2" = list(
        Treatment = "none",
        Age = 68.1,
        Sex = "M",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Pembro, Ipi",
        Time_from_last_Tx = 1,
        Class = "On ICI Therapy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = "alive",
        Case = "Control",
        Class_short = "ICI",
        Other_irAEs = "None",
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    ),
    "NC3" = list(
        Treatment = "Ipilimumab",
        Age = 70.5,
        Sex = "F",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Pembro, Durva",
        Time_from_last_Tx = 12,
        Class = "On ICI Therapy",
        Diagnosis = "enteritis",
        CTCAE = 2,
        MES = 0,
        RECIST_response = "PD",
        Overall_survival = 15.1,
        Status = "alive",
        Case = "Control",
        Class_short = "ICI",
        Other_irAEs = "irEnteritis",
        Drug = "Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4)",
        DrugType = "Monotherapy"
    ),
    "NC4" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 63.7,
        Sex = "M",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "none",
        Time_from_last_Tx = NA,
        Class = "On ICI Therapy",
        Diagnosis = "normal",
        CTCAE = 1,
        MES = 0,
        RECIST_response = "PR",
        Overall_survival = 12.6,
        Status = "alive",
        Case = "Control",
        Class_short = "ICI",
        Other_irAEs = "None",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "NC5" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 58.5,
        Sex = "M",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Nivo, Ipi+Nivo",
        Time_from_last_Tx = 1,
        Class = "On ICI Therapy",
        Diagnosis = "enteritis",
        CTCAE = 1,
        MES = 0,
        RECIST_response = "PD",
        Overall_survival = 12.3,
        Status = "deceased",
        Case = "Control",
        Class_short = "ICI",
        Other_irAEs = "irEnteritis",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "NC6" = list(
        Treatment = "Ipilimumab + Nivolumab",
        Age = 49.3,
        Sex = "F",
        Malignancy = "Melanoma",
        Stage = 4,
        Prior_Tx = "Melanoma",
        Time_from_last_Tx = 2,
        Class = "On ICI Therapy",
        Diagnosis = "normal",
        CTCAE = 2,
        MES = 0,
        RECIST_response = "PD",
        Overall_survival = 10.4,
        Status = "alive",
        Case = "Control",
        Class_short = "ICI",
        Other_irAEs = "None",
        Drug = "Anti-PD-1, Anti-CTLA-4",
        DrugName = "Ipilimumab (Anti-CTLA4),Nivolumab (anti-PD1)",
        DrugType = "Combination"
    ),
    "CT1" = list(
        Treatment = NA,
        Age = 21.4,
        Sex = "F",
        Malignancy = "melanoma",
        Stage = "1A",
        Prior_Tx = NA,
        Time_from_last_Tx = NA,
        Class = "Healthy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = NA,
        Case = "Control",
        Class_short = "H",
        Other_irAEs = NA,
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    ),
    "CT2" = list(
        Treatment = NA,
        Age = 61.7,
        Sex = "M",
        Malignancy = NA,
        Stage = NA,
        Prior_Tx = NA,
        Time_from_last_Tx = NA,
        Class = "Healthy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = NA,
        Case = "Control",
        Class_short = "H",
        Other_irAEs = NA,
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    ),
    "CT3" = list(
        Treatment = NA,
        Age = 50.5,
        Sex = "M",
        Malignancy = NA,
        Stage = NA,
        Prior_Tx = NA,
        Time_from_last_Tx = NA,
        Class = "Healthy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = NA,
        Case = "Control",
        Class_short = "H",
        Other_irAEs = NA,
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    ),
    "CT4" = list(
        Treatment = NA,
        Age = 43.4,
        Sex = "M",
        Malignancy = NA,
        Stage = NA,
        Prior_Tx = NA,
        Time_from_last_Tx = NA,
        Class = "Healthy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = NA,
        Case = "Control",
        Class_short = "H",
        Other_irAEs = NA,
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    ),
    "CT5" = list(
        Treatment = NA,
        Age = 70.6,
        Sex = "F",
        Malignancy = "gastric",
        Stage = 2,
        Prior_Tx = NA,
        Time_from_last_Tx = NA,
        Class = "Healthy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = NA,
        Case = "Control",
        Class_short = "H",
        Other_irAEs = NA,
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    ),
    "CT6" = list(
        Treatment = NA,
        Age = 44.3,
        Sex = "M",
        Malignancy = NA,
        Stage = NA,
        Prior_Tx = NA,
        Time_from_last_Tx = NA,
        Class = "Healthy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = NA,
        Case = "Control",
        Class_short = "H",
        Other_irAEs = NA,
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    ),
    "CT7" = list(
        Treatment = NA,
        Age = 50.2,
        Sex = "M",
        Malignancy = NA,
        Stage = NA,
        Prior_Tx = NA,
        Time_from_last_Tx = NA,
        Class = "Healthy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = NA,
        Case = "Control",
        Class_short = "H",
        Other_irAEs = NA,
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    ),
    "CT8" = list(
        Treatment = NA,
        Age = 61.5,
        Sex = "F",
        Malignancy = NA,
        Stage = NA,
        Prior_Tx = NA,
        Time_from_last_Tx = NA,
        Class = "Healthy",
        Diagnosis = "normal",
        CTCAE = 0,
        MES = 0,
        RECIST_response = NA,
        Overall_survival = NA,
        Status = NA,
        Case = "Control",
        Class_short = "H",
        Other_irAEs = NA,
        Drug = "None",
        DrugName = NA,
        DrugType = NA
    )
)



cd_clusters <- list(
    "CD3.cell.annotation.txt" = c("1" = "MAIT", "2" = "Type 3 cytokines Trm", "3" = "Type 1 cytokines Trm", "4" = "Naive/CM", "5" = "5", "6" = "T follicular helper", "7" = "Treg", "8" = "Th1 effector", "9" = "Cycling", "10" = "Cytotoxic Effector", "11" = "LP Trm", "12" = "CD8 Trm", "13" = "IEL: CD8/gdT"),
    "CD4.cell.annotation.txt" = c("1" = "Th17 Trm", "2" = "Th1 Trm-A", "3" = "Th1 Trm-R", "4" = "Naive", "5" = "Central Memory", "6" = "Tfh", "7" = "Treg", "8" = "Cytotoxic", "9" = "Th1 Effector", "10" = "Cycling"),
    "CD8.cell.annotation.txt" = c("1" = "Trm - IEL", "2" = "Trm - LP1", "3" = "Trm - LP2", "4" = "MAIT", "5" = "Term effector", "6" = "CM/Naive", "7" = "Cytotoxic Effector", "8" = "Cycling")
)

cd_clusters[["CD3.cell.annotation.txt"]]
