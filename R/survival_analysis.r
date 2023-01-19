# Purpose
# Survival analysis of clusters

# Takes in a dataframe with survival data and returns a survival curve fitted
# against variables specified in factor_colname.

#' Fit survival curves by a grouping factor
#'
#' This function takes a dataframe that contains survival data, specifically
#' columns for survival time and survival censor, as well as a grouping
#' variable that is detailed in a column, factor_colname.
#'
#' @param survival_df A dataframe that contains survival information
#' @param survival_time_colname A string literal that refers to a column
#' in survival_df that outlines patient survival time
#' @param survival_censor_colname A string literal that refers to a column
#' in survial_df that outlines patient survival censoring
#' @param factor_colname A string literal that refers to a column in
#' survival_df that contains a grouping variable used to categorize patients.
#' Survival between these groups will be measured.
#'
#' @return An object with fitted survival curves
#' @export
fit_survival_curve_by_factor <- function(
    survival_df,
    survival_time_colname,
    survival_censor_colname,
    factor_colname
    ) {
        # Formatting survival formula
        surv_formula <- stats::as.formula(
            base::paste0(
                "Surv(",
                survival_time_colname,
                ",",
                survival_censor_colname,
                ") ~ ",
                factor_colname
            )
        )

        fit <- survminer::surv_fit(
            surv_formula,
            data = survival_df
        )

        return(fit)
}

#' Plot survival curves
#'
#' This function takes a set of fitted survival curves and plots them.
#'
#' @param survival_fit An object of fitted survival curves
#' @param survival_df A dataframe containing survival data used to fit
#' the survival curves in survival_fit
#' @param factor_colname A string literal that refers to a column in
#' survival_df that contains grouping variables to categorize patients.
#'
#' @return A kaplan-meier survival plot.
#' @export
plot_survplot <- function(survival_fit, survival_df, factor_colname) {
    survminer::ggsurvplot(
        data = survival_df,
        fit = survival_fit,
        pval = TRUE,
        palette = "jco",
        xlab = "Months",
        ylab = "Survival probability [%]"
    )
}

#' Find median cutoff
#'
#' Adds a column that determines if a sample is above or below median
#' in a given metric
#'
#' @param input_df The dataframe containing the metric
#' @param metric_colname A string literal that refers to a column in input_df
#'
#' @return input_df but with an extra column detailing median status
#' @export
assign_median_cutoff <- function(
    input_df,
    metric_colname
    ) {
        # Find median
        median_stat <- stats::median(input_df[[metric_colname]])

        # Adding cutoff
        result <- input_df %>%
            dplyr::mutate(
                median_status = ifelse(
                    !!dplyr::sym(metric_colname) > median_stat,
                    "above",
                    "below"
                )
            )
        return(result)
    }

#' Find maxstat survival cutoff
#'
#' Adds a column that determines if a metric is above or below the
#' maxstat cutoff.
#'
#' @param input_df The dataframe containing the metric
#' @param metric_colname A string literal that refers to a column in input_df
#' @param censor_colname A string literal referring to a column of censors
#' @param survival_time_colname A string literal referring to a column for survival time
#'
#' @return input_df but with an extra column detailing median status
#' @export
assign_maxstat_cutoff <- function(
    input_df,
    metric_colname,
    censor_colname,
    survival_time_colname
    ) {
        # Finding maxstat cutoff
        form <- paste0(
                    survival_time_colname,
                    ",",
                    censor_colname
                )

        # Assigning maxstat
        maxstat_cutoff <- maxstat::maxstat.test(
            survival::Surv(
                paste0(form, ") ~", !!dplyr::sym(metric_colname)),
                data = input_df,
                smethod = "LogRank",
                pmethod = "Lau94"
            )
        )

        # Assigning cutoff to input_df
        result <- input_df %>%
            dplyr::mutate(
                maxstat_status = ifelse(
                    !!dplyr::sym(metric_colname) > maxstat_cutoff,
                    "above",
                    "below"
                )
            )
        return(result)
    }