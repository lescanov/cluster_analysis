# Purpose:
# Perform differential gene expression analysis between two groups

library(edgeR)
library(dplyr)
library(tibble)
library(stats)
library(conflicted)
library(org.Hs.eg.db)

# org.Hs.eg.db overwrites select from dplyr
conflict_prefer_all("dplyr")

# This function creates an array of factors which have each factor repeated
# by the number of times they are present in factor_column.
# The inputs ref_factor and compare_factor are elements of factor_column.
create_normalization_factors <- function(
    df_with_factor,
    factor_colname = "cluster",
    ref_factor,
    compare_factor) {
        df <- df_with_factor %>%
            filter(!!sym(factor_colname) %in% c(ref_factor, compare_factor)) %>%
            group_by(!!sym(factor_colname)) %>%
            summarise(n = n())

        length_ref <- df %>%
            filter(!!sym(factor_colname) %in% ref_factor) %>%
            select(n) %>%
            deframe()

        length_compare <- df %>%
            filter(!!sym(factor_colname) %in% compare_factor) %>%
            select(n) %>%
            deframe()

        normalization_factors <- c(
            factor(rep(ref_factor, length_ref)),
            factor(rep(compare_factor, length_compare))
        )

        normalization_factors <- relevel(
            normalization_factors,
            ref = ref_factor
        )

        return(normalization_factors)
}

# Extracts patient IDs from df_with_factor and then arrange raw_counts by
# patient IDs of the corresponding factors.
arrange_df_by_norm_factors <- function(
    raw_counts,
    df_with_factor,
    ref_factor,
    compare_factor,
    factor_colname = "cluster",
    id_colname = "patient_id") {
        # Geting patient IDs for each factor
        ref_ids <- df_with_factor %>%
            filter(!!sym(factor_colname) %in% ref_factor) %>%
            select(all_of(id_colname)) %>%
            deframe()

        compare_ids <- df_with_factor %>%
            filter(!!sym(factor_colname) %in% compare_factor) %>%
            select(all_of(id_colname)) %>%
            deframe()

        # Arranging counts based on order of factors (referene -> compare)
        arranged_counts <- raw_counts %>%
            select(all_of(c(ref_ids, compare_ids)))

        return(arranged_counts)
}

# Identifies differentially expressed genes between two groups.
# input_df is a dataframe which contains raw counts where patients are
# ordered by their respective normalization factor.
# Output is a dataframe containing differentially expressed genes
# filtered by FDR values < 0.01.
find_degs <- function(input_df, normalization_factors) {
    # Convert to matrix for edgeR
    input_matrix <- as.matrix(input_df)

    # Normalizing with TMM
    counts <- DGEList(input_matrix)
    keep_counts <- rowSums(counts$counts) > 50
    counts <- counts[keep_counts, , keep.lib.sizes = FALSE]
    counts <- calcNormFactors(counts, method = "TMM")

    # Creating model matrix using normalization factors
    model_matrix <- model.matrix(~normalization_factors)

    # Find DEGs using glmQLF tests
    estimated_dispersion <- estimateDisp(counts, model_matrix, robust = TRUE)
    fit <- glmQLFit(estimated_dispersion, model_matrix)
    ftest <- glmQLFTest(fit)
    result <- topTags(ftest, n = nrow(input_df))$table

    # Filter DEGs for significant results
    result <- result %>%
        filter(FDR <0.01) %>%
        arrange(desc(logFC)) %>%
        rownames_to_column(var = "ensembl_id")

    return(result)
}

# Converts gene symbols into 
ensembl_to_entrez_and_symbol <- function(input_df, id_colname = "ensembl_id") {
    # Converting ensembl to gene symbol
    gene_ids <- mapIds(
        org.Hs.eg.db,
        keys = input_df[[id_colname]],
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
    )

    # Converting ensembl to entrez ids
    entrez_ids <- mapIds(
        org.Hs.eg.db,
        keys = input_df[[id_colname]],
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
    )

    # Appending converted ids
    result <- input_df %>%
        mutate(
            gene_symbol = gene_ids,
            entrez_id = entrez_ids
        )

    return(result)
}

find_degs_between_clusters <- function(
    df_with_factor,
    factor_colname = "cluster",
    ref_factor,
    compare_factor,
    raw_counts,
    id_colname = "patient_id"
) {
    norm_factors <- create_normalization_factors(
        df_with_factor = df_with_factor,
        factor_colname = factor_colname,
        ref_factor = ref_factor,
        compare_factor = compare_factor
    )

    arranged_counts <- arrange_df_by_norm_factors(
        raw_counts = raw_counts,
        df_with_factor = df_with_factor,
        ref_factor = ref_factor,
        compare_factor = compare_factor,
        factor_colname = factor_colname,
        id_colname = id_colname
    )

    degs <- find_degs(arranged_counts, norm_factors) %>%
        ensembl_to_entrez_and_symbol()

    return(degs)
}
