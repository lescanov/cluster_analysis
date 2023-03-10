% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_expression_analysis.r
\name{create_normalization_factors}
\alias{create_normalization_factors}
\title{Create normalization factors for edgeR}
\usage{
create_normalization_factors(
  df_with_factor,
  factor_colname = "cluster",
  ref_factor,
  compare_factor
)
}
\arguments{
\item{df_with_factor}{Dataframe with a column corresponding to
categorical factors}

\item{factor_colname}{A string literal that species the column with
categorical factors}

\item{ref_factor}{The factor that will be used a reference in the
differential expression analysis}

\item{compare_factor}{The factor that will be compared against
the reference factor}
}
\value{
Returns a factor with levels corresponding to ref_factor and
compare_factor, which reflect the number of samples assigned to each factor.
}
\description{
Create normalization factors for edgeR
}
\examples{
data(iris)
norm_factors <- create_normalization_factors(
     iris,
     factor_colname = "Species",
     ref_factor = "setosa",
     compare_factor = "virginica"
)
}
