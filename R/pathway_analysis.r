# Purpose
# Perform pathway analysis, specifically GSEA and ORA.

library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(ggplot2)

# Perform GSEA analysis on a sorted vector of logFCs, named with entrez IDs.
# Can perform GSEA using GO (using  All modules), KEGG or Reactome.
# df_of_degs is intended to be the output of find_degs_between_clusters,
# and this function assumes there are columns called logFC and entrez_id
# in the input dataframe.
# pathway annotation must be one of go, kegg or reactome.
# output is a dataframe listing pathways and pvalues.
perform_gsea <- function(df_of_degs, pathway_annotation) {
    if (!pathway_annotation %in% c("go", "kegg", "reactome")) {
        return(print("insert one of go, kegg or reactome"))
    }

    # Removing any NAs in df_of_degs
    df <- df_of_degs %>%
        filter(is.na(entrez_id) == FALSE)

    # For GSEA, must supply a sorted vector of logFCs with names as entrez ids
    gene_list <- df[["logFC"]]
    names(gene_list) <- df[["entrez_id"]]

    if (pathway_annotation == "go") {
        result <- gseGO(
            geneList = gene_list,
            OrgDb = "org.Hs.eg.db",
            ont = "ALL",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH"
        )
    } else
    
    if (pathway_annotation == "kegg") {
        result <- gseKEGG(
            geneList = gene_list,
            keyType = "ncbi-geneid",
            pvalueCutoff = 0.05
        )
    } else

    if (pathway_annotation == "reactome") {
        result <- gsePathway(
            geneList = gene_list
        )
    }

    return(result)
}

# Perform ORA analysis on a vector of entrez IDs.
# Can perform ORA using GO (using  All modules), KEGG or Reactome.
# This specifically will perform pathway analysis on up or downregulated
# genes, which must be specified by the up_or_down parameter.
# df_of_degs is intended to be the output of find_degs_between_clusters,
# and this function assumes there are columns called logFC and entrez_id
# in the input dataframe.
# pathway annotation must be one of go, kegg or reactome.
# output is a dataframe listing pathways and pvalues.
perform_ora_on_deregulated_genes <- function(
    df_of_degs,
    pathway_annotation,
    up_or_down
    ) {
    if (!pathway_annotation %in% c("go", "kegg", "reactome")) {
        return(print("insert one of go, kegg or reactome"))
    }

    if (!up_or_down %in% c("up", "down")) {
        return(print("insert one of up or down"))
    }

    # Removing any NAs in df_of_degs
    df <- df_of_degs %>%
        filter(is.na(entrez_id) == FALSE)

    # Extracting up/downregulated genes
    if (up_or_down == "up") {
        upregulated_genes <- df %>%
            filter(logFC > 0)
        gene_list <- upregulated_genes[["entrez_id"]]
    } else {
        downregulated_genes <- df %>%
            filter(logFC < 0)
        gene_list <- downregulated_genes[["entrez_id"]]
    }

    # Performing pathway analysis on extracted genes
    if (pathway_annotation == "go") {
        result <- enrichGO(
            gene = gene_list,
            ont = "ALL",
            OrgDb = org.Hs.eg.db
      )
    } else

    if (pathway_annotation == "kegg") {
        result <- enrichKEGG(
            gene = gene_list,
            keyType = "ncbi-geneid",
            qvalueCutoff = 0.05
        )
    } else

    if (pathway_annotation == "reactome") {
        result <- enrichPathway(
            gene = gene_list,
            qvalueCutoff = 0.05
        )
    }

    return(result)
}

# Save dotplot (shows first 10 results) automatically to results folder
save_dotplot <- function(pathway_analysis, title) {
    plot <- dotplot(pathway_analysis, title =  title)
    ggsave(
        plot = plot,
        filename = paste0(title, ".pdf"),
        path = "results/",
        height = 6,
        width = 8,
        units = "in"
    )
}

# Save ridgeplot for GSEA analyses (does not work for ORA), automatically
# to the results folder.
save_ridgeplot <- function(pathway_analysis, title) {
    plot <- ridgeplot(pathway_analysis, showCategory = 10)
    ggsave(
        plot = plot,
        filename = paste0(title, ".pdf"),
        path = "results/",
        height = 6,
        width = 8,
        units = "in"
    )
}