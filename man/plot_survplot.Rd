% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_analysis.r
\name{plot_survplot}
\alias{plot_survplot}
\title{Plot survival curves}
\usage{
plot_survplot(survival_fit, survival_df, factor_colname)
}
\arguments{
\item{survival_fit}{An object of fitted survival curves}

\item{survival_df}{A dataframe containing survival data used to fit
the survival curves in survival_fit}

\item{factor_colname}{A string literal that refers to a column in
survival_df that contains grouping variables to categorize patients.}
}
\value{
A kaplan-meier survival plot.
}
\description{
This function takes a set of fitted survival curves and plots them.
}
