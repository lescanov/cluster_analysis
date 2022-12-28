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

#' Create multiple dataframes with subset vs whole(without subset)
#'
#' This creates a list of multiple different dataframes, consisting of
#' subset vs whole (without subset) columns.
#'
#' @param input_df A dataframe that contains: a column with a grouping variable
#' from which one can extract a subset of the data and a column containing
#' a metric that the user wishes to compare between subset and rest of dataset.
#' @param grouping_column The column that contains the grouping variable, such
#' as patient subtype, etc.
#' @param metric_column The column that contains a metric the user would like
#' to compare between two subsets of the dataset.
#' @param list_of_subset_groups A list of string literals specifying subgroups
#' of the grouping column.
#'
#' @return A list of adataframes of the format of define_subset_against_whole()
#' @export
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

#' Verify equality of variances in samples grouped by a variable
#'
#' This function uses Levene's test to determine if variance is equal
#' amongst all subgroups.
#'
#' @param input_df A dataframe that contains: a column with a grouping variable
#' from which one can extract a subset of the data and a column containing
#' a metric that the user wishes to compare between subset and rest of dataset.
#' @param grouping_column The column that contains the grouping variable, such
#' as patient subtype, etc.
#' @param metric_column The column that contains a metric the user would like
#' to compare between two subsets of the dataset.
#'
#' @return A boolean value, TRUE if the p-value obtained from Levene's test is
#' greater than 0.05 and FALSE if the p-value is less tha 0.05.
#' @export
verify_equality_of_variances <- function(
    input_df,
    grouping_colname,
    metric_colname
) {
    test_formula <- stats::as.formula(
        paste0(metric_colname, "~", grouping_colname)
    )

    # Testing ewuality of variances between the two groups
    levene_result <- rstatix::levene_test(input_df, test_formula)[["p"]]

    if (levene_result > 0.05) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Verify if distribution amongst grouping variables is normal
#'
#' This function determines if the distribution of a particular metric
#' is normal, across grouping variables, using a Shapiro-Wilk test. It
#' is expected that this function will be used to assess which test is
#' valid for the subset vs whole (without subset) comparisons.
#'
#' @param input_df A dataframe that contains: a column with a grouping variable
#' from which one can extract a subset of the data and a column containing
#' a metric that the user wishes to compare between subset and rest of dataset.
#' @param grouping_column The column that contains the grouping variable, such
#' as patient subtype, etc.
#' @param metric_column The column that contains a metric the user would like
#' to compare between two subsets of the dataset.
#' @param subset A string literal that specifies a subgroup in the grouping
#' column
#'
#' @return A boolean value, TRUE if the p-value obtained from Shapiro-Wilk's
#' test is greater than 0.05 and FALSE if the p-value is less tha 0.05.
#' @export
verify_normality <- function(
    input_df,
    grouping_colname,
    metric_colname,
    subset = "subset"
) {
    # Getting subsets of data
    subset_data <- input_df %>%
        dplyr::filter(!!dplyr::sym(grouping_colname) %in% subset)

    without_subset_data <- input_df %>%
        dplyr::filter(!(!!dplyr::sym(grouping_colname)) %in% subset)

    # Testing that both subsets of dataset have normal distribution
    shapiro_subset <- subset_data %>%
        rstatix::shapiro_test(!!dplyr::sym(metric_colname))
    shapiro_without_subset <- without_subset_data %>%
        rstatix::shapiro_test(!!dplyr::sym(metric_colname))

    if (shapiro_subset > 0.05 && shapiro_without_subset > 0.05) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Verify if distributions under different grouping variables are equal
#'
#' This function is meant to test equality of distributions, to determine
#' if the Mann-Whitney test is appropriate. It uses the Kolomogorov-Smirnov test
#' to assess if two distributions are equal.
#'
#' @param input_df A dataframe that contains: a column with a grouping variable
#' from which one can extract a subset of the data and a column containing
#' a metric that the user wishes to compare between subset and rest of dataset.
#' @param grouping_column The column that contains the grouping variable, such
#' as patient subtype, etc.
#' @param metric_column The column that contains a metric the user would like
#' to compare between two subsets of the dataset.
#' @param subset A string literal that specifies a subgroup in the grouping
#' column
#'
#' @return A boolean value, TRUE if the p-value obtained from Kolmogorov-Smirnov
#' test is greater than 0.05 and FALSE if the p-value is less tha 0.05.
#' @export
verify_equal_disitributions <- function(
    input_df,
    grouping_colname,
    metric_colname,
    subset = "subset"
) {
    # Defining data for subset
    subset_data <- input_df %>%
        dplyr::filter(!!dplyr::sym(grouping_colname) %in% subset) %>%
        dplyr::select(!!dplyr::sym(metric_colname)) %>%
        tibble::deframe()

    without_subset_data <- input_df %>%
        dplyr::filter(!(!!dplyr::sym(grouping_colname)) %in% subset) %>%
        dplyr::select(!!dplyr::sym(metric_colname)) %>%
        tibble::deframe()

    # Perform Kolmogorov-Smirnov test
    ks_result <- ks.test(
        x = without_subset_data,
        y = subset_data,
        data = input_df
    )

    if (ks_result > 0.05) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

susbet_vs_whole_t_test <- function(
    input_df,
    grouping_colname = "group",
    metric_colanme,
    equal_variance = FALSE
) {
    # Check if equal_variance is a boolean value
    stopifnot(is.logical(equal_variance))

    # Create comparison formula using metric and grouping variables
    comparison_formula <- stats::as.formula(
        paste0(metric_colanme, "~", grouping_colname)
    )

    # Perform t-test, depending on value of equal variance, will perform
    # etiher Welch's (equal variance = FALSE) or Students (TRUE)
    result <- input_df %>%
        rstatix::t_test(
            comparison_formula,
            var.equal = equal_variance
        )

    return(result)
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

adjusted_subset_vs_whole_test <- function(
    input_df,
    grouping_colname,
    metric_colname,
    list_of_subset_groups,
    test_to_use
) {
    stopifnot(test_to_use %in% c("student", "welch", "mann_whitney"))

    # Create a list of dataframes that utilize multiple different subsets
    # of the same dataset
    list_of_subset_comparisons <- create_list_of_subset_comparisons(
        input_df = input_df,
        grouping_colname = grouping_colname,
        metric_colname = metric_colname,
        list_of_subset_groups = list_of_subset_groups
    )

    if (test_to_use == "mann_whitney") {
        test_result <- purrr::map(
            list_of_subset_comparisons,
            subset_vs_whole_mann_whitney,
            metric_colname = metric_colname
        )
    } else if (test_to_use == "welch") {
        test_result <- purrr::map(
            list_of_subset_comparisons,
            subset_vs_whole_t_test,
            metric_colname = metric_colname,
            equal_variance = FALSE
        )
    } else if (test_to_use == "student") {
        test_result <- purrr::map(
            list_of_subset_comparisons,
            subset_vs_whole_t_test,
            metric_colname = metric_colname,
            equal_variance = TRUE
        )
    }

    # Combining list into one dataframe then correcting p-value
    result_df <- test_result %>%
        dplyr::bind_rows() %>%
        rstatix::adjust_pvalue(method = "fdr")

    return(result_df)
}

#' Determine which hypothesis test is appropriate for distributions
#'
#' There are three hypothesis tests that are considered for this function:
#' Student t-test (assumption of variance is equal), Welch's t-test
#' (assumotion of unequal variances) when distribution of residuals
#' are normal. When normality is not achieved, this function
#' will assess if the Mann-Whitney U-test is appropriate by testing
#' if distributions are equal via. the Kolomogrov-Smirnov test.
#'
#' @param input_df A dataframe that contains: a column with a grouping variable
#' from which one can extract a subset of the data and a column containing
#' a metric that the user wishes to compare between subset and rest of dataset.
#' @param grouping_column The column that contains the grouping variable, such
#' as patient subtype, etc.
#' @param metric_column The column that contains a metric the user would like
#' to compare between two subsets of the dataset.
#' @param subset A string literal that specifies a subgroup in the grouping
#' column
#'
#' @return A string literal of: "mann_whitney", "welch", "student" or "other".
#' These represent the appropriate hypothesis test to use given the
#' distributions of the subset and whole (without subset) portions of the
#' dataset. The "other" represents a case where all the tests fail and
#' an alternative test should be used.
#' @export
decide_which_hypothesis_test <- function(
    input_df,
    grouping_colname,
    metric_colname,
    subset
) {
    normal_distribution <- verify_normality(
        input_df = input_df,
        grouping_colname = grouping_colname,
        metric_colname = metric_colname
    )

    if (normal_distribution == FALSE) {
        print("Shapiro-Wilk test failed: Non-normal distribution")
        equal_distribution <- verify_equal_disitributions(
            input_df = input_df,
            grouping_colname = grouping_colname,
            metric_colname = metric_colname
        )

        if (equal_distribution == TRUE) {
            result <- "mann_whitney"
            print("Kolmogorov-Smnirov test passed: Mann-Whitney is appropriate")
            return(result)
        } else {
            result <- "other"
            print("Kolmogorov-Smnirov test failed: Must find alternative test")
            return(result)
        }
    }
    print("Shapiro-Wilk's test passed: Normal distribution assumed")

    equal_variance <- verify_equality_of_variances(
        input_df = input_df,
        grouping_colname = grouping_colname,
        metric_colanme = metric_colname
    )

    if (equal_variance == FALSE) {
        result <- "welch"
        print("Levene's test failed: Use Welch's t-test")
    } else {
        result <- "student"
        print("Leven's test passed: Use student's t-test")
    }
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