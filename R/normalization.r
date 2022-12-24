# Purpose:
# Normalize raw-seq reads based on selected subtype

# Retrieve patient ids from a dataframe with a specified subtype.
# subtype_to_filter is an element found in the subtype_column.
# The output is an array of patient_ids of class character

#' Select Patient/Sample IDs based on subtype
#'
#' @param input_df A dataframe which contains a column specifying
#' patient subtypes
#' @param id_column The column where sample/patient IDs are stored
#' @param subtype_column The column where patient subtypes are stored
#' @param subtype_to_filter A stirng literal represting a subtype
#' present in subtype_column.
#'
#' @return A vector of characters specifying IDs for a respective subtype
#' @export
select_ids_from_subtype <- function(
    input_df, id_column, subtype_column, subtype_to_filter) {
    ids <- input_df %>%
        dplyr::filter({{subtype_column}} %in% subtype_to_filter) %>%
        dplyr::select(dplyr::all_of({{id_column}})) %>%
        tibble::deframe() %>%
        base::as.character()

    return(ids)
}

#' Normalize raw counts using edgeR
#'
#' This function takes a dataframe of raw RNA-seq counts (samples as columns),
#' and normalizes it using edgeR's TMM method.
#' It returns a dataframe of log2 transofrmed normalized counts.
#' This function will select a subset of patients.
#' This functions offers the option of liberal or strict
#' count filtering. Liberal filtering sees if counts per gene
#' are greater than 50. The strict method utilizes edgeR's
#' filterByExpr function and asserts that a cpm of 10 (default value)
#' is achieved in a significant amount of samples.
#'
#' @param input_df A dataframe of raw RNA-seq counts, samples are columns.
#' @param ids_to_use Patient IDs used to select a subset of patients.
#' @param filter_method A string literal of either "strict" or "liberal",
#' liberal will take genes with rowSums > 50, strict will use edgeR's
#' filterByExpr which takes a cpm of 10 (by default) in a significant
#' number of samples.
#'
#' @return A dataframe of log2 transformed TMM normalized counts,
#' where patients are columns, genes are rows
#' @export
normalize_counts <- function(input_df, ids_to_use, filter_method = "liberal") {
    stopifnot(filter_method %in% c("strict", "liberal"))

    # Select patients of subtype to normalize
    df <- input_df %>%
        dplyr::select(dplyr::all_of(ids_to_use))

    # Conver to matrix for edgeR
    input_matrix <- base::as.matrix(df)
    counts <- edgeR::DGEList(input_matrix)

    # Applying filter method to counts
    if (filter_method == "strict") {
        keep_counts <- edgeR::filterByExpr(counts)
    } else
    if (filter_method == "liberal") {
        keep_counts <- stats::rowSums(counts$counts) > 50
    }

    # Normalizing with TMM
    counts <- counts[keep_counts, , keep.lib.sizes = FALSE]
    counts <- edgeR::calcNormFactors(counts, method = "TMM")
    counts <- edgeR::cpm(counts, log = TRUE)
    counts <- counts %>%
        base::as.data.frame() %>%
        tibble::rownames_to_column(var = "ensembl_id")

    return(counts)
}

# Assuming using output of normalize_counts, this will extract the expression
# of a gene and pivot the df to have columns gene, patient_id and expression

#' Extract gene expression from normalize_counts() result
#'
#' This function takes a dataframe which is the result from normalize_counts(),
#' and pivots it such that there are now only two columns, patient_ids and
#' log2_expression. The input dataframe is filtered for gene(s) of interest,
#' and then pivoted.
#'
#' @param normalized_counts A dataframe from normalize_counts() where
#' samples are columns, genes are rwos.
#' @param gene_column A column in normalized_counts that specifies where
#' the gene symbols are stored.
#' @param gene A string literal representing the gene(s) of interest.
#'
#' @return A dataframe with two columns, patient_ids and log2_expression
#' @export
extract_gene_expression <- function(normalized_counts, gene_column, gene) {
    df <- normalized_counts %>%
        dplyr::filter({{gene_column}} %in% gene) %>%
        tidyr::pivot_longer(
            !{{gene_column}},
            names_to = "patient_id",
            values_to = "log2_expression"
        )

    return(df)
}
