# For analyzing expression across user defined subtypes.
# These tools are meant for measuring gene expression across cytogenetic
# and molecular subtypes.
# Kruskal wallis is for cytogenetic subtypes
# A modified version of the Mann-whitney U test is used for molecular
# subtypes in order to maintain assumption of independence between comparisons.
# This takes a subset of the total sample population (e.g. mutation)
# and compares this to the total population without the subset.

# October 31, 2022

library(rstatix)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(purrr)

# Performs and plots kruskal-wallis. Performs Dunn's test as post-hoc test.
# p-values are adjusted for multiple comparisons. Plots boxplots with
# the x_variable as the numeric variable and the y_variable as the groupping variable.
# Coordinates are flipped for formatting.
perform_and_plot_kruskal_test <- function(df, expression, group) {
    # Defining formula with column names
    form <- as.formula(paste0(expression, "~", group))

    # Performing kruskal wallis
    kruskal_result <- df %>%
        kruskal_test(form)

    # Performing dunn test's as post hoc test
    dunn_test <- df %>%
        dunn_test(
            form,
            p.adjust.method = "fdr"
        ) %>%
        add_xy_position(
            x = group
        )

    # Plotting with boxplot
    boxplot <- ggboxplot(
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
    stat_pvalue_manual(dunn_test, hide.ns = TRUE) +
    labs(
        subtitle = get_test_label(kruskal_result, detailed = TRUE),
        caption = get_pwc_label(dunn_test)
    ) +
    theme(legend.position = "none")

    print(boxplot)
}

# The following tools are for using the modified mann-whitney.
# The workflow of the following tools is to work on lists of subsets
# and necessitates that each function is mapped to each element of the list.
# Every variable with the "_column" suffix is expected to be a non-string
# name of the column in input_df.
# identifier_column is the column containing the different subtypes
# sample_column is the sample IDs to the sequencing samples and is needed
# to ensure that each patient is included only once in the test sets.
define_subset_comparisons <- function(input_df, subset_group, reference_group, identifier_column, sample_column) {
    # Defining the subset group
    subset_df <- input_df %>%
        filter({{identifier_column}} %in% subset_group) %>%
        distinct({{sample_column}}, .keep_all = TRUE)

    # Identify patients that are different among reference and subset
    subset_patients <- subset_df %>%
        select({{sample_column}}) %>%
        deframe()

    # Defining reference group
    reference_df <- input_df %>%
        filter({{identifier_column}} %in% reference_group) %>%
        distinct({{sample_column}}, .keep_all = TRUE)

    # Defining reference, without subset patients
    reference_without_subset <- reference_df %>%
        filter(!{{sample_column}} %in% subset_patients)

    # Binding reference_without_subset and subset
    reference_and_subset <- reference_without_subset %>%
        bind_rows(subset_df) %>%
        mutate(group = paste(subset_group, "vs", reference_group))

    return(reference_and_subset)
}

perform_subset_wilcox <- function(input_df, identifier_colname, gene) {
    comparison_formula <- as.formula(paste0(gene, "~", identifier_colname))
    result <- input_df %>%
        wilcox_test(comparison_formula, paired = FALSE, alternative = "two.sided")
    return(result)
}

plot_subset_comparisons <- function(comparisons_df, x_variable, y_variable) {
    ggboxplot(
        comparisons_df,
        x_variable,
        y_variable,
        palette = "jco",
        add = "jitter",
        color = x_variable,
        repel = TRUE,
        xlab = ""
    ) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme(legend.position = "none")

}

# This function utillzes a list of dataframes specifying expression of a y_variable
# across different subtypes specified by x_varbiable.
# Using the define_subset_comparisons function will produce a column called group
# that details which comparison is being made.
plot_each_subset_comparison <- function(x_variable = "subtype", y_variable, list_of_subtypes) {
    # Creating a list of plots
    list_of_plots <- map(
        list_of_subtypes,
        plot_subset_comparisons,
        x_variable = x_variable,
        y_variable = y_variable
    )

    # Combining plots into one plot
    plot <- ggarrange(plotlist = list_of_plots)

    return(plot)
}