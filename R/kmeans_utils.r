# Purpose:
# Utility functions for downstream analysis of kmeans clustering

library(dplyr)
library(ggplot2)
library(factoextra)
library(tibble)
library(ggpubr)

# Utilize three tests: elbow, silhouette and gap to determine optimal
# amount of clusters to use. Prints all three plots
plot_optimal_cluster_tests <- function(input_df) {
    print(fviz_nbclust(input_df, kmeans, "wss"))
    print(fviz_nbclust(input_df, kmeans, "gap_stat"))
    print(fviz_nbclust(input_df, kmeans, "silhouette"))
}

# Extracts clusters from kmeans_result object
extract_clusters_from_kmeans <- function(kmeans_cluster_result) {
    # Extract clusters from kmeans cluster result
    clusters <- as.data.frame(kmeans_cluster_result[["cluster"]]) %>%
        rownames_to_column(var = "patient_id") %>%
        rename("cluster" = contains("cluster"))

    # Convert clusters to factor
    clusters[["cluster"]] <- as.factor(clusters[["cluster"]])

    return(clusters)
}

# Performs kmeans and then extracts column
perform_kmeans_and_extract_clusters <- function(
    random_seed,
    input_df,
    centers,
    nstart = 25
) {
    # Seed random seed
    set.seed(random_seed)

    # Perform kmeans
    kmeans_result <- kmeans(input_df, centers = centers, nstart = nstart)

    # Printing cluster results
    print(fviz_cluster(kmeans_result, data = input_df))

    result <- extract_clusters_from_kmeans(kmeans_result)

    return(result)
}