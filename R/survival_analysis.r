# Purpose
# Survival analysis of clusters

library(survival)
library(survminer)

# Takes in a dataframe with survival data and returns a survival curve fitted
# against variables specified in factor_colname.
fit_survival_curve_by_factor <- function(
    survival_df,
    survival_time_colname,
    survival_censor_colname,
    factor_colname
    ) {
        # Formatting survival formula
        surv_formula <- as.formula(
            paste0(
                "Surv(",
                survival_time_colname,
                ",",
                survival_censor_colname,
                ") ~ ",
                factor_colname
            )
        )

        fit <- surv_fit(
            surv_formula,
            data = survival_df
        )

        return(fit)
}

plot_survplot <- function(survival_fit, survival_df, factor_colname) {
    ggsurvplot(
        data = survival_df,
        fit = survival_fit,
        pval = TRUE,
        palette = "jco",
        xlab = "Months",
        ylab = "Survival probability [%]"
    )
}