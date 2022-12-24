# Purpose
# Perform PCA on raw counts and extract the first 2 principal components

#' Perform Principal Component analysis on normalized RNA-seq counts
#'
#' @param normalized_counts Normalized RNA-seq counts, where columns are samples
#' @param gene_id_colname The column name in normalized_counts for gene ids
#'
#' @return A prcomp object of centred and scaled principal components, with patients as rows and genes as columns
#' @export
#'
#' @examples
#' rna_seq_counts <- some_data_frame
#' pca_result <- perform_pca_on_counts(rna_seq_counts)
perform_pca_on_counts <- function(normalized_counts, gene_id_colname) {
    df <- normalized_counts %>%
        tibble::column_to_rownames(var = gene_id_colname) %>%
        t()

    # Perform PCA
    pca <- stats::prcomp(df, scale = TRUE, center = TRUE)

    return(pca)
}

#' Extract specified amount of principal components from a prcomp object
#'
#' @param prcomp_result A prcomp result, designed to be result from perform_pca_on_counts
#' @param number_of_pc The number of principal components to be extracted from prcomp_result
#'
#' @return A dataframe containing the specified amount of principal components, with PCs as columns
#' @export
#'
#' @examples
#' # Extract first 2 principal components from prcomp result
#' pca <- perform_pca_on_counts(some_rnaseq_counts)
#' pc_1_and_2 <- extract_principal_components(pca, number_of_pc = 2)
extract_principal_components <- function(prcomp_result, number_of_pc) {
    # Principal components are stored as value x in prcomp object
    # Convert to dataframe in order to select components
    principal_components <- as.data.frame(prcomp_result[["x"]])

    # Select number of components
    principal_components <- principal_components %>%
        dplyr::select(dplyr::all_of(1:number_of_pc))

    return(principal_components)
}

# Performs PCA on counts, prints plot of variance explained by principal
# components, and extracts principal components from resulting prcomp result.

#' Perform PCA and extract specified amount of principal components
#'
#' @param normalized_counts Normalized RNA-seq counts, where columns are samples
#' @param gene_id_colname The column name in normalized_counts for gene ids
#' @param number_of_pc The number of principal components to be extracted from prcomp_result
#'
#' @return A dataframe containing the specified amount of principal components, with PCs as columns
#' @export
#'
#' @examples
#' # Extract first 2 principal components from normalized counts
#' pc_1_and_2 <- perform_and_extract_pca(some_ranseq_counts, gene_id_colname, 2)
perform_and_extract_pca <- function(
    normalized_counts,
    gene_id_colname,
    number_of_pc_to_extract
) {
    # Performing CPA
    pca_result <- perform_pca_on_counts(
        normalized_counts = normalized_counts,
        gene_id_colname = gene_id_colname
    )

    # Print variance
    print(factoextra::fviz_eig(pca_result))

    # Extract clusters
    result <- extract_principal_components(
        prcomp_result = pca_result,
        number_of_pc = number_of_pc_to_extract
    )

    return(result)
}

#' Inner join principal components to a dataframe
#'
#' @param df_to_append A dataframe to append principal components
#' @param principal_component_df A dataframe with principal components
#' @param colname_to_append_by A string specifying a column that is present in both df_to_append and principal_component_df
#'
#' @return A single dataframe that has principal components appended to df_to_append
#' @export
#'
#' @examples
#' pc_1_2 <- perform_and_extract_pca(rnaseq_counts, gene_id, 2)
#' df <- some_data_frame
#' df_with_pc_1_2 <- append_principal_components(
#'      df,
#'      pc_1_2,
#'      "some_colname"
#' )
append_principal_components <- function(
    df_to_append,
    principal_component_df,
    colname_to_append_by) {
        stopifnot(
            nrow(df_to_append) == nrow(principal_component_df),
            colname_to_append_by %in% colnames(df_to_append) & colname_to_append_by %in% colnames(principal_component_df)
        )

        # Extract components
        pc_df <- principal_component_df %>%
            tibble::rownames_to_column(var = colname_to_append_by)

        # Append to input_df
        result <- df_to_append %>%
            dplyr::inner_join(
                pc_df,
                by = colname_to_append_by
            )

        return(result)
}
