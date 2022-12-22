# Purpose:
# Normalize raw-seq reads based on selected subtype

library(edgeR)
library(dplyr)
library(tibble)
library(tidyr)

# Retrieve patient ids from a dataframe with a specified subtype.
# subtype_to_filter is an element found in the subtype_column.
# The output is an array of patient_ids of class character
select_ids_from_subtype <- function(
    input_df, id_column, subtype_column, subtype_to_filter) {
    ids <- input_df %>%
        filter({{subtype_column}} %in% subtype_to_filter) %>%
        select({{id_column}}) %>%
        deframe() %>%
        as.character()

    return(ids)
}

# Normalize raw reads using edgeR (TMM method)
# Have two different filter methods, strict and liberal.
# Strict uses the function filterByExpr for filtering whilst
# liberal uses rowSums > 50.
#
# input_df format must be patients as columns, genes as rownames
# output is dataframe with patients as columns, with a column for gene ids
# Counts are log2 transformed.
#
# ids_to_use are patient ids that correspond to subtype that
# will be normalized.
normalize_counts <- function(input_df, ids_to_use, filter_method = "strict") {
    if (!filter_method %in% c("strict", "liberal")) {
        return(print("insert strict or liberal for filter method"))
    }

    # Select patients of subtype to normalize
    df <- input_df %>%
        select(all_of(ids_to_use))

    # Conver to matrix for edgeR
    input_matrix <- as.matrix(df)
    counts <- DGEList(input_matrix)

    # Applying filter method to counts
    if (filter_method == "strict") {
        keep_counts <- filterByExpr(counts)
    } else
    if (filter_method == "liberal") {
        keep_counts <- rowSums(counts$counts) > 50
    }

    # Normalizing with TMM
    counts <- counts[keep_counts, , keep.lib.sizes = FALSE]
    counts <- calcNormFactors(counts, method = "TMM")
    counts <- cpm(counts, log = TRUE)
    counts <- counts %>%
        as.data.frame() %>%
        rownames_to_column(var = "ensembl_id")

    return(counts)
}

# Assuming using output of normalize_counts, this will extract the expression
# of a gene and pivot the df to have columns gene, patient_id and expression
extract_gene_expression <- function(normalized_counts, gene_column, gene) {
    df <- normalized_counts %>%
        filter({{gene_column}} %in% gene) %>%
        pivot_longer(
            !{{gene_column}},
            names_to = "patient_id",
            values_to = "log2_expression"
        )

    return(df)
}
