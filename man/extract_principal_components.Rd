% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.r
\name{extract_principal_components}
\alias{extract_principal_components}
\title{Extract specified amount of principal components from a prcomp object}
\usage{
extract_principal_components(prcomp_result, number_of_pc)
}
\arguments{
\item{prcomp_result}{A prcomp result, designed to be result from perform_pca_on_counts}

\item{number_of_pc}{The number of principal components to be extracted from prcomp_result}
}
\value{
A dataframe containing the specified amount of principal components, with PCs as columns
}
\description{
Extract specified amount of principal components from a prcomp object
}
\examples{
# Extract first 2 principal components from prcomp result
# Do not run
# pca <- perform_pca_on_counts(some_rnaseq_counts)
# pc_1_and_2 <- extract_principal_components(pca, number_of_pc = 2)
}
