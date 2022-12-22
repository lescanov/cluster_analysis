# Purpose
# Perform PCA on raw counts and extract the first 2 principal components

library(dplyr)
library(tibble)

# Performs PCA (scaled and centered) on normalized seq counts.
# Expected input is dataframe of counts where columns are patients and
# have gene_ids in a column specified by gene_id_colname.
# The function converts gene_id column to
# rownames then transposes the dataframe (where patients are now rows,
#  genes as columns) then performs PCA.
perform_pca_on_counts <- function(normalized_counts, gene_id_colname) {
    df <- normalized_counts %>%
        column_to_rownames(var = gene_id_colname) %>%
        t()

    # Perform PCA
    pca <- prcomp(df, scale = TRUE, center = TRUE)

    return(pca)
}

# Extract specified amount of principal components from a prcomp_result object.
# The result is a dataframe with principal components as columns and
# patient_ids as rownames.
extract_principal_components <- function(prcomp_result, number_of_pc) {
    # Principal components are stored as value x in prcomp object
    # Convert to dataframe in order to select components
    principal_components <- as.data.frame(prcomp_result[["x"]])

    # Select number of components
    principal_components <- principal_components %>%
        select(all_of(1:number_of_pc))

    return(principal_components)
}

# Performs PCA on counts, prints plot of variance explained by principal
# components, and extracts principal components from resulting prcomp result.
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
    print(fviz_eig(pca_result))

    # Extract clusters
    result <- extract_principal_components(
        prcomp_result = pca_result,
        number_of_pc = number_of_pc_to_extract
    )

    return(result)
}

# Inner join specified amount of principal components to a dataframe.
append_principal_components <- function(
    df_to_append,
    principal_component_df,
    colname_to_append_by) {
        # Extract components
        pc_df <- principal_component_df %>%
            rownames_to_column(var = colname_to_append_by)

        # Append to input_df
        result <- df_to_append %>%
            inner_join(
                pc_df,
                by = colname_to_append_by
            )

        return(result)
}
