% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_expression_analysis.r
\name{find_degs}
\alias{find_degs}
\title{Identify differentially expression genes bewteen two groups}
\usage{
find_degs(arranged_raw_counts, normalization_factors)
}
\arguments{
\item{arranged_raw_counts}{A dataframe of raw counts, with columns as
samples that arranged by group, and genes as rownames. It is assumed that
gene IDs will be in ENSEMBL format.}

\item{normalization_factors}{A factor with two levels, reference and
comparison.}
}
\value{
A dataframe containing differentially expressed genes
between the two groups. Genes that have a FDR p=value of less than
0.01 are retained, all others are filtered out.
The test used is Genewise Negative Binomial Generalized Linear Models
with Quasi-Likelihood tests, from the glmQLFIT test.
The resulting dataframe also has a column reserved for the ENSEMBL IDs.
}
\description{
Identify differentially expression genes bewteen two groups
}
