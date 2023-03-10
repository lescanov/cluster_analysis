% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathway_analysis.r
\name{perform_ora_on_deregulated_genes}
\alias{perform_ora_on_deregulated_genes}
\title{Perform Over-Representation Analysis on differentially expressed genes}
\usage{
perform_ora_on_deregulated_genes(df_of_degs, pathway_annotation, up_or_down)
}
\arguments{
\item{df_of_degs}{A dataframe of differentially expressed genes.
Must contain a column for entrez geneIDs and logFC. Anticipated to
be a result of edgeR's topTags funcion.}

\item{pathway_annotation}{Must be a string literal matching one of
"kegg", "go" or "reactome". Specifies which pathway annotation to
use for pathway analysis.}

\item{up_or_down}{A string literal of either "up" or "down", specifying
if the user would like to perform ORA on up or downregulated genes.}
}
\value{
An enrichResult object
}
\description{
This function utilizes clusterProfiler's algorithim for ORA. Specifically,
this will perform ora on a subset of differentially regulated genes,
either upregulated or downregulated genes. The purpose of this design
choice is to compliment GSEA's analysis by analyzing which
particular pathways are being upregulated or downregulated specifically.
The user can choose 3 different pathway annotation databases:
"go", "kegg" or "reactome". For Gene Ontology, all databases are used.
}
