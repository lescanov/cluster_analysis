library(dplyr)
library(readr)
library(readxl)
library(stringr)

# My criteria for general classification of cytogenetic subgroups:
# > 1. For non-complex karyotype, aberration must be seen in all clones
# > 2. For complex karyotype, majority of clones had complex karyotypes
# > 3. Complex karyotypes contained 3+ genomic lesions
# > 4. For normal karyotypes, had consulted the othterCytogenetics column
#      If the karyotype was normal, but FISH had detected a translocation,
#      did not consider patient as normal karyotype.

# Currently have only cytogenetics assigned for patients used in manuscript
# Need to expand my list to be more comprehensive and include all patients
# First identifying which patients do not have cytogenetics assigned
# Importing patient data
clinical_summary <- read_excel("data/beataml_wv1to4_clinical.xlsx", sheet = 1)
manuscript_cytogenetics <- read_csv("data/cytogenetics_assigned.csv")

# First assign cytogenetics based on karyotype
partial_cytogenetics <- clinical_summary %>%
  filter(rnaSeq %in% "y") %>%
  mutate(
    cytogenetic_subtype = case_when(
      str_detect(karyotype, "^46,X[XY]\\[..\\]") ~ "Normal karyotype",
      str_detect(consensusAMLFusions, "PML-RARA") ~ "t(15;17)",
      str_detect(consensusAMLFusions, "MLLT3-KMT2A|KMT2A_re") ~ "MLL rearranged", #nolint
      str_detect(consensusAMLFusions, "CBFB-MYH11") ~ "inv(16)",
      str_detect(consensusAMLFusions, "RUNX1-RUNX1T1") ~ "t(8;21)"
    ),
    vitalStatus = ifelse(vitalStatus == "Dead", 1, 0)
  )

patients_with_assigned_cytogenetics <- partial_cytogenetics %>%
    filter(is.na(cytogenetic_subtype) == FALSE) %>%
    select(
        dbgap_subject_id,
        karyotype,
        cytogenetic_subtype,
        otherCytogenetics
    )

# Identifying which patients that:
# 1) Have RNAseq
# 2) Do not have cytogenetics assigned
no_cytogenetics_assigned <- partial_cytogenetics %>%
    filter(!dbgap_subject_id %in% manuscript_cytogenetics[["dbgap_subject_id"]]) %>% #nolint
    filter(is.na(cytogenetic_subtype) == TRUE) %>%
    select(
        dbgap_subject_id,
        karyotype,
        cytogenetic_subtype,
        otherCytogenetics
    )

write.csv(
    no_cytogenetics_assigned,
    file = "data/patients_with_no_cytogenetic_subtype.csv"
)

# Manually assigned rest of the patients with no cytogenetics
assigned_cytogenetics <- read_csv("data/non_manuscript_cytogenetics_assigned.csv") #nolint

# Merging patient cytogenetics
# This is currently missing patients found in partial cytogenetics
manually_assigned_cytogenetics <- manuscript_cytogenetics %>%
    bind_rows(assigned_cytogenetics) %>%
    bind_rows(patients_with_assigned_cytogenetics)

# Appendining manually assigned cytogenetics with rest of clinical summary
complete_clinical_summary <- partial_cytogenetics %>%
    filter(dbgap_subject_id %in% manually_assigned_cytogenetics[["dbgap_subject_id"]]) %>%
    inner_join(manually_assigned_cytogenetics, by = "dbgap_subject_id")

colnames(complete_clinical_summary)[97:101] <- gsub(
    "\\..*",
    "",
    colnames(complete_clinical_summary)[97:101]
)

write.csv(complete_clinical_summary, "data/clinical_summary.csv")
