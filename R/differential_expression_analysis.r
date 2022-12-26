# Purpose:
# Perform differential gene expression analysis between two groups

#' Create normalization factors for edgeR
#'
#' @param df_with_factor Dataframe with a column corresponding to
#' categorical factors
#' @param factor_colname A string literal that species the column with
#' categorical factors
#' @param ref_factor The factor that will be used a reference in the
#' differential expression analysis
#' @param compare_factor The factor that will be compared against
#' the reference factor
#'
#' @return Returns a factor with levels corresponding to ref_factor and
#' compare_factor, which reflect the number of samples assigned to each factor.
#' @export
#'
#' @examples
#' data(iris)
#' norm_factors <- create_normalization_factors(
#'      iris,
#'      factor_colname = "Species",
#'      ref_factor = "setosa",
#'      compare_factor = "virginica"
#' )
create_normalization_factors <- function(
    df_with_factor,
    factor_colname = "cluster",
    ref_factor,
    compare_factor) {
        df <- df_with_factor %>%
            dplyr::filter(!!dplyr::sym(factor_colname) %in% c(ref_factor, compare_factor)) %>%
            dplyr::group_by(!!dplyr::sym(factor_colname)) %>%
            dplyr::summarise(n = dplyr::n())

        length_ref <- df %>%
            dplyr::filter(!!dplyr::sym(factor_colname) %in% ref_factor) %>%
            dplyr::select(n) %>%
            tibble::deframe()

        length_compare <- df %>%
            dplyr::filter(!!dplyr::sym(factor_colname) %in% compare_factor) %>%
            dplyr::select(n) %>%
            tibble::deframe()

        normalization_factors <- c(
            base::factor(rep(ref_factor, length_ref)),
            base::factor(rep(compare_factor, length_compare))
        )

        normalization_factors <- stats::relevel(
            normalization_factors,
            ref = ref_factor
        )

        return(normalization_factors)
}

#' Arrange raw counts by the order of specified factors for differential expression analysis
#' 
#' @param raw_counts A dataframe of raw counts for use in expression analysis,
#' with sample/patient IDs as columns
#' @param df_with_factor A dataframe that contains a column of factors to
#' identify groups for expression analysis
#' @param ref_factor A string literal for the reference factor found in factor
#' @param compare_factor A stirng literal for the factor contrasted against
#' colname the reference factor, found in factor_colname
#' @param factor_colname A string literal specifying a column in df_with_factor
#' @param id_colname A string literal specifying a column that contains
#' sample/patient IDs
#'
#' @return A dataframe of raw counts with columns ordered by a specified set of factors
#' @export
arrange_df_by_norm_factors <- function(
    raw_counts,
    df_with_factor,
    ref_factor,
    compare_factor,
    factor_colname = "cluster",
    id_colname = "patient_id") {
        # Geting patient IDs for each factor
        ref_ids <- df_with_factor %>%
            dplyr::filter(!!dplyr::sym(factor_colname) %in% ref_factor) %>%
            dplyr::select(dplyr::all_of(id_colname)) %>%
            tibble::deframe()

        compare_ids <- df_with_factor %>%
            dplyr::filter(!!dplyr::sym(factor_colname) %in% compare_factor) %>%
            dplyr::select(dplyr::all_of(id_colname)) %>%
            tibble::deframe()

        # Arranging counts based on order of factors (referene -> compare)
        arranged_counts <- raw_counts %>%
            dplyr::select(dplyr::all_of(c(ref_ids, compare_ids)))

        return(arranged_counts)
}

#' Identify differentially expression genes bewteen two groups
#'
#' @param arranged_raw_counts A dataframe of raw counts, with columns as
#' samples that arranged by group, and genes as rownames. It is assumed that
#' gene IDs will be in ENSEMBL format.
#' @param normalization_factors A factor with two levels, reference and
#' comparison.
#'
#' @return A dataframe containing differentially expressed genes
#' between the two groups. Genes that have a FDR p=value of less than
#' 0.01 are retained, all others are filtered out.
#' The test used is Genewise Negative Binomial Generalized Linear Models
#' with Quasi-Likelihood tests, from the glmQLFIT test.
#' The resulting dataframe also has a column reserved for the ENSEMBL IDs.
#' @export
find_degs <- function(arranged_raw_counts, normalization_factors) {
    # Convert to matrix for edgeR
    input_matrix <- base::as.matrix(arranged_raw_counts)

    # Normalizing with TMM
    counts <- edgeR::DGEList(input_matrix)
    keep_counts <- base::rowSums(counts$counts) > 50
    counts <- counts[keep_counts, , keep.lib.sizes = FALSE]
    counts <- edgeR::calcNormFactors(counts, method = "TMM")

    # Creating model matrix using normalization factors
    model_matrix <- stats::model.matrix(~normalization_factors)

    # Find DEGs using glmQLF tests
    estimated_dispersion <- edgeR::estimateDisp(counts, model_matrix, robust = TRUE)
    fit <- edgeR::glmQLFit(estimated_dispersion, model_matrix)
    ftest <- edgeR::glmQLFTest(fit)
    result <- edgeR::topTags(ftest, n = nrow(arranged_raw_counts))$table

    # Filter DEGs for significant results
    result <- result %>%
        dplyr::filter(FDR <0.01) %>%
        dplyr::arrange(dplyr::desc(logFC)) %>%
        tibble::rownames_to_column(var = "ensembl_id")

    return(result)
}

#' Convert ENSEMBL IDs into gene symbols and entrez-IDs
#'
#' @param deg_df A dataframe that is a result of the topTags edgeR function.
#' Can be a dataframe returned from find_degs function.
#' @param id_colname A string literal that specifies the column in the dataframe
#' reserved for ENSEMBL IDs. By default, assumed to be "ensembl_id".
#'
#' @return deg_df with two new columns, one for gene symbols, another for
#' entrez IDs
#' @export
ensembl_to_entrez_and_symbol <- function(deg_df, id_colname = "ensembl_id") {
    # Converting ensembl to gene symbol
    gene_ids <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = deg_df[[id_colname]],
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
    )

    # Converting ensembl to entrez ids
    entrez_ids <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = deg_df[[id_colname]],
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
    )

    # Appending converted ids
    result <- deg_df %>%
        dplyr::mutate(
            gene_symbol = gene_ids,
            entrez_id = entrez_ids
        )

    return(result)
}

#' Identify differentially expressed genes between two groups, while generating normalization factors, arranged raw counts, and columns of converted IDs.
#'
#' @param raw_counts A dataframe of raw counts for use in expression
#' analysis, with sample/patient IDs as columns
#' @param df_with_factor A dataframe that contains a column of factors
#' to identify groups for expression analysis
#' @param ref_factor A string literal for the reference factor found in
#' factor_colname
#' @param compare_factor A stirng literal for the factor contrasted against
#' the reference factor, found in factor_colname
#' @param factor_colname A string literal specifying a column in df_with_factor
#' @param id_colname A string literal specifying a column that contains
#' sample/patient IDs
#'
#' @return A dataframe resulting from the edgeR topTags function, detailing
#' differentially expressed genes between two groups, with columns for gene
#' symbols and entrez IDs.
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
