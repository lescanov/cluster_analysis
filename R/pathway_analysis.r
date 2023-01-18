# Purpose
# Perform pathway analysis, specifically GSEA and ORA.

#' Perform Gene Set Enrichment Analysis on differentially expressed genes
#'
#' Utilize clusterProfiler's GSEA algorithim to perform pathway analysis
#' on a dataframe of differentially expressed genes. It is expected
#' that the input dataframe of differentially expressed genes has
#' a column for entrez gene IDs and a column for log fold change called
#' "logFC".
#' The user can choose 3 different pathway annotation databases:
#' KEGG, Gene Ontology and Reactome. For Gene Ontology, all the
#' databases are used.
#'
#' @param df_of_degs A dataframe of differentially expressed genes.
#' Must contain a column for entrez geneIDs and logFC. Anticipated to
#' be a result of edgeR's topTags funcion.
#' @param pathway_annotation Must be a string literal matching one of
#' "kegg", "go" or "reactome". Specifies which pathway annotation to
#' use for pathway analysis.
#'
#' @return A gseaResult object
#' @export
perform_gsea <- function(df_of_degs, pathway_annotation) {
    base::stopifnot(
        pathway_annotation %in% c("go", "kegg", "reactome"),
        "entrez_id" %in% colnames(df_of_degs)
    )

    # Removing any NAs in df_of_degs
    df <- df_of_degs %>%
        dplyr::filter(base::is.na(entrez_id) == FALSE)

    # For GSEA, must supply a sorted vector of logFCs with names as entrez ids
    gene_list <- df[["logFC"]]
    names(gene_list) <- df[["entrez_id"]]

    if (pathway_annotation == "go") {
        result <- clusterProfiler::gseGO(
            geneList = gene_list,
            OrgDb = "org.Hs.eg.db",
            ont = "ALL",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH"
        )
    } else

    if (pathway_annotation == "kegg") {
        result <- clusterProfiler::gseKEGG(
            geneList = gene_list,
            keyType = "ncbi-geneid",
            pvalueCutoff = 0.05
        )
    } else

    if (pathway_annotation == "reactome") {
        result <- ReactomePA::gsePathway(
            geneList = gene_list
        )
    }

    return(result)
}

#' Perform Over-Representation Analysis on differentially expressed genes
#'
#' This function utilizes clusterProfiler's algorithim for ORA. Specifically,
#' this will perform ora on a subset of differentially regulated genes,
#' either upregulated or downregulated genes. The purpose of this design
#' choice is to compliment GSEA's analysis by analyzing which
#' particular pathways are being upregulated or downregulated specifically.
#' The user can choose 3 different pathway annotation databases:
#' "go", "kegg" or "reactome". For Gene Ontology, all databases are used.
#'
#' @param df_of_degs A dataframe of differentially expressed genes.
#' Must contain a column for entrez geneIDs and logFC. Anticipated to
#' be a result of edgeR's topTags funcion.
#' @param pathway_annotation Must be a string literal matching one of
#' "kegg", "go" or "reactome". Specifies which pathway annotation to
#' use for pathway analysis.
#' @param up_or_down A string literal of either "up" or "down", specifying
#' if the user would like to perform ORA on up or downregulated genes.
#'
#' @return An enrichResult object
#' @export
perform_ora_on_deregulated_genes <- function(
    df_of_degs,
    pathway_annotation,
    up_or_down
    ) {
    base::stopifnot(
        pathway_annotation %in% c("go", "kegg", "reactome"),
        up_or_down %in% c("up", "down"),
        "entrez_id" %in% colnames(df_of_degs)
    )

    # Removing any NAs in df_of_degs
    df <- df_of_degs %>%
        dplyr::filter(base::is.na(entrez_id) == FALSE)

    # Extracting up/downregulated genes
    if (up_or_down == "up") {
        upregulated_genes <- df %>%
            dplyr::filter(logFC > 0)
        gene_list <- upregulated_genes[["entrez_id"]]
    } else {
        downregulated_genes <- df %>%
            dplyr::filter(logFC < 0)
        gene_list <- downregulated_genes[["entrez_id"]]
    }

    # Performing pathway analysis on extracted genes
    if (pathway_annotation == "go") {
        result <- clusterProfiler::enrichGO(
            gene = gene_list,
            ont = "ALL",
            OrgDb = org.Hs.eg.db
      )
    } else

    if (pathway_annotation == "kegg") {
        result <- clusterProfiler::enrichKEGG(
            gene = gene_list,
            keyType = "ncbi-geneid",
            qvalueCutoff = 0.05
        )
    } else

    if (pathway_annotation == "reactome") {
        result <- ReactomePA::enrichPathway(
            gene = gene_list,
            qvalueCutoff = 0.05
        )
    }

    return(result)
}

# Save dotplot (shows first 10 results) automatically to results folder

#' Automatically save dotplot of pathway analysis
#'
#' @param pathway_analysis A gseaResult or enrichResult object
#' @param title Title of the file to be saved
#' @param path Path of the file destination. Defaults to "results/"
#'
#' @return A 8in (width) x 6in (height) pdf file. Dotplot shows top 10 pathways.
#' @export
save_dotplot <- function(pathway_analysis, title, path = "results/") {
    plot <- clusterProfiler::dotplot(pathway_analysis, title =  title)
    ggplot2::ggsave(
        plot = plot,
        filename = paste0(title, ".pdf"),
        path = path,
        height = 6,
        width = 8,
        units = "in"
    )
}

perform_gsea_and_ora <- function()

# Save ridgeplot for GSEA analyses (does not work for ORA), automatically
# to the results folder.

#' Automatically save ridgeplot of pathway analysis
#'
#' NOTE: This function does not work with ORA analyses.
#'
#' @param pathway_analysis A gseaResult object
#' @param title Title of the file to be saved
#' @param path Path of the file destination. Defaults to "results/"
#' @param number_of_pathways An unsigned number that specifies how many
#' pathways to show in the ridgeplot. Defaults to 15.
#'
#' @return A ridgeplot of 8in (width) x 6in (height) pdf file. Ridgeplot shows
#' specified number of pathways. Defaults to 15.
#' @export
save_ridgeplot <- function(
    pathway_analysis,
    title,
    path = "results/",
    number_of_pathways = 15
    ) {
    plot <- clusterProfiler::ridgeplot(
        pathway_analysis,
        showCategory = number_of_pathways
    )
    ggplot2::ggsave(
        plot = plot,
        filename = paste0(title, ".pdf"),
        path = path,
        height = 6,
        width = 8,
        units = "in"
    )
}