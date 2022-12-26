# For analyzing expression across user defined subtypes.
# These tools are meant for measuring gene expression across cytogenetic
# and molecular subtypes.
# Kruskal wallis is for cytogenetic subtypes
# A modified version of the Mann-whitney U test is used for molecular
# subtypes in order to maintain assumption of independence between comparisons.
# This takes a subset of the total sample population (e.g. mutation)
# and compares this to the total population without the subset.

#' Perform and plot a Kruskal-Wallis test
#'
#' Given a grouping variable, this function will perform a Kruskal-Wallis
#' test, along with Dunn's test as a post-hoc test.
#'
#' @param df A dataframe that contains the grouping variable and the
#' metric to compare across groups.
#' @param expression A string literal for the column representing gene
#' expression.
#' @param group A string literal for the column representing the
#' grouping variable.
#'
#' @return A boxplot comparing expression across groups. Plots significance
#' baed on results of Dunn's test. The p-value obtained from the Kruskal Wallis
#' along with number of samples compared and degrees of freedom will be reported
#' at the top of the plot.
#' @export
perform_and_plot_kruskal_test <- function(df, expression, group) {
    # Defining formula with column names
    form <- stats::as.formula(paste0(expression, "~", group))

    # Performing kruskal wallis
    kruskal_result <- df %>%
        rstatix::kruskal_test(form)

    # Performing dunn test's as post hoc test
    dunn_test <- df %>%
        rstatix::dunn_test(
            form,
            p.adjust.method = "fdr"
        ) %>%
        rstatix::add_xy_position(
            x = group
        )

    # Plotting with boxplot
    boxplot <- ggpubr::ggboxplot(
        df,
        x = group,
        y = expression,
        color = group,
        palette = "jco",
        add = "jitter",
        size = 1.2,
        xlab = "",
        repel = TRUE
    ) +
    ggpubr::stat_pvalue_manual(dunn_test, hide.ns = TRUE) +
    ggplot2::labs(
        subtitle = rstatix::get_test_label(kruskal_result, detailed = TRUE),
        caption = rstatix::get_pwc_label(dunn_test)
    ) +
    ggplot2::theme(legend.position = "none")

    print(boxplot)
}

#' Define groups: Subset of dataset and Dataset without subset
#'
#' This function separates samples into two groups. The first, a subset of
#' samples (such as a group of patients belonging to a particular substype)
#' and the second, the rest of the dataset without samples of this subset.
#' This function assumes that there is a column specifying sample groupings,
#' and that there is both a grouping for the subset the user is interested in
#' as well as a grouping for all samples in the dataset.
#' This function is used to set hypotheses test downstream, if the user
#' wishes to compare a subset against an "All" group. This is particularly
#' useful especially when samples can inhabit multiple different subtypes
#' for example a cancer patient can harbour multiple mutations. If
#' the user is interested to see if the expression of a gene is enriched
#' in a particular subset of patients, this can serve as a statistical
#' basis in testing if this is actually so.
#' The main purppose of this function is maintain the assumption of
#' independence that many inferrential hypothesis tests require.
#'
#' @param input_df A dataframe that contains samples the user wishes to compare.
#' This dataframe should have a grouping for both the subset and the whole
#' dataset.
#' @param subset_group A string literal that defines a subset of samples the
#' user wishes to compare to the rest of the dataset.
#' @param reference_group A string literal that defines samples belonging to the
#' whole dataset.
#' @param identifier_column A column that contains the grouping variables
#' subset_group and reference_group.
#' @param sample_column A column that contains the sample IDs for this dataset.
#'
#' @return A dataframe with a column that specifies a comparison between
#' the outlined subset and the rest of the dataset.
#' @export
#'
#' @examples
#' # Example for data preparation prior to using the function
#' data(iris)
#' # Defining a grouping that encapsulates the entire dataset
#' all_group <- iris %>% dplyr::mutate(group = "all")
#'
#' # Defining a subset of the dataset
#' sestosa_group <- iris %>%
#'      dplyr::filter(Species %in% "sestosa") %>%
#'      dplyr::mutate(group = "sestosa")
#'
#' # Combining dataframes, this can then be used in the function.
#' # The group column will be identifier_column.
#' whole_and_subset <- rbind(all_group, sestosa_group)
define_subset_comparisons <- function(input_df, subset_group, reference_group, identifier_column, sample_column) {
    # Defining the subset group
    subset_df <- input_df %>%
        dplyr::filter({{identifier_column}} %in% subset_group) %>%
        dplyr::distinct({{sample_column}}, .keep_all = TRUE)

    # Identify patients that are different among reference and subset
    subset_patients <- subset_df %>%
        dplyr::select({{sample_column}}) %>%
        tibble::deframe()

    # Defining reference group
    reference_df <- input_df %>%
        dplyr::filter({{identifier_column}} %in% reference_group) %>%
        dplyr::distinct({{sample_column}}, .keep_all = TRUE)

    # Defining reference, without subset patients
    reference_without_subset <- reference_df %>%
        dplyr::filter(!{{sample_column}} %in% subset_patients)

    # Binding reference_without_subset and subset
    reference_and_subset <- reference_without_subset %>%
        dplyr::bind_rows(subset_df) %>%
        dplyr::mutate(group = paste(subset_group, "vs", reference_group))

    return(reference_and_subset)
}

#' Create a dataframe with grouping variables subset vs whole (without subset)
#'
#' The purpose of this function is to create a dataframe where one wants to
#' compare a subset of a dataset against the rest of the dataset. The main
#' reason to do this is to maintain the assumption of indepdence in
#' many inferrential hypohtesis tests. This would allow the user
#' to see how much a metric changes in the subset, relative to the
#' rest of the data in the dataset.
#'
#' @param input_df A dataframe that contains: a column with a grouping variable
#' from which one can extract a subset of the data and a column containing
#' a metric that the user wishes to compare between subset and rest of dataset.
#' @param grouping_column The column that contains the grouping variable, such
#' as patient subtype, etc.
#' @param subset_group A string literal present in grouping_column that defines
#' the subset subgroup the user is interested in.
#' @param metric_column The column that contains a metric the user would like
#' to compare between two subsets of the dataset.
#'
#' @return A dataframe with 3 columns: grouping_column, metric_column and group,
#' which is a column that defines which part of the data belongs to the subset
#' and the rest of the dataset.
#' @export
define_subset_against_whole <- function(
    input_df,
    grouping_colname,
    subset_group,
    metric_colname
    ) {
    # Select just grouping and metric column
    df <- input_df %>%
        dplyr::select(
            !!dplyr::sym(grouping_colname),
            !!dplyr::sym(metric_colname)
            )

    # Select subset
    subset <- df %>%
        dplyr::filter(!!dplyr::sym(grouping_colname) %in% subset_group) %>%
        dplyr::mutate(group = "subset")

    # Define the dataset without the subset subtype
    whole_without_subset <- df %>%
        dplyr::filter(!(!!dplyr::sym(grouping_colname)) %in% subset_group) %>%
        dplyr::mutate(group = "whole_without_subset")

    # Combining these two together into one dataframe
    subset_against_whole <- dplyr::bind_rows(subset, whole_without_subset)

    return(subset_against_whole)
}

create_list_of_subset_comparisons <- function(
    input_df,
    grouping_colname,
    metric_colname,
    list_of_subset_groups
) {
    # Create a list of dataframes that utilize multiple different subsets
    # of the same dataset
    list_of_subset_comparisons <- purrr::map(
        list_of_subset_groups,
        define_subset_against_whole,
        input_df = input_df,
        grouping_column = grouping_colname,
        metric_column = metric_colname
    )

    return(list_of_subset_comparisons)
}

#' Compare a subset against the rest of dataset using Mann-Whitney U Tests.
subset_vs_whole_mann_whitney <- function(
    input_df,
    grouping_colname = "group",
    metric_colname
    ) {
    # Creating a formula to compare metric across groups
    comparison_formula <- stats::as.formula(
        paste0(metric_colname, "~", grouping_colname)
    )

    # Performing Mann-Whitney U-test
    result <- input_df %>%
        rstatix::wilcox_test(
            comparison_formula, paired = FALSE, alternative = "two.sided"
        )

    return(result)
}

adjusted_subset_vs_whole_mann_whitney <- function(
    input_df,
    grouping_colname,
    metric_colname,
    list_of_subset_groups
) {
    # Create a list of dataframes that utilize multiple different subsets
    # of the same dataset
    list_of_subset_comparisons <- create_list_of_subset_comparisons(
        input_df = input_df,
        grouping_colname = grouping_colname,
        metric_colname = metric_colname,
        list_of_subset_groups = list_of_subset_groups
    )

    # Supplying this list to Mann whitney
    list_of_mann_whitney_results <- purrr::map(
        list_of_subset_comparisons,
        subset_vs_whole_mann_whitney,
        grouping_colname = grouping_colname,
        metric_colname = metric_colname
    )

    # Combining list into one dataframe then correcting p-value
    result_df <- list_of_mann_whitney_results %>%
        dplyr::bind_rows() %>%
        rstatix::adjust_pvalue(method = "fdr")

    return(result_df)
}


plot_subset_comparisons <- function(
    input_df,
    grouping_colname,
    metric_colname
    ) {
    ggpubr::ggboxplot(
        input_df,
        grouping_colname,
        metric_colname,
        palette = "jco",
        add = "jitter",
        color = grouping_colname,
        repel = TRUE,
        xlab = ""
    ) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::theme(legend.position = "none")

}

# This function utillzes a list of dataframes specifying expression of a y_variable
# across different subtypes specified by x_varbiable.
# Using the define_subset_comparisons function will produce a column called group
# that details which comparison is being made.
plot_each_subset_comparison <- function(
    grouping_colname,
    metric_colname,
    list_of_subset_groups,
    input_df
    ) {
    list_of_subset_comparisons <- create_list_of_subset_comparisons(
        input_df = input_df,
        grouping_colname = grouping_colname,
        metric_colname = metric_colname,
        list_of_subset_groups = list_of_subset_groups
    )

    # Creating a list of plots
    list_of_plots <- purrr::map(
        list_of_subset_comparisons,
        plot_subset_comparisons,
        grouping_colname = grouping_colname,
        metric_colname = metric_colname
    )

    # Combining plots into one plot
    plot <- ggpubr::ggarrange(plotlist = list_of_plots)

    return(plot)
}