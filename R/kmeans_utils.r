# Purpose:
# Utility functions for downstream analysis of kmeans clustering

# Utilize three tests: elbow, silhouette and gap to determine optimal
# amount of clusters to use. Prints all three plots

#' Plot elbow, gap statistic and silhouette plots to determine optimal amount
#' of clusters for kmeans clustering
#'
#' This function will print plots for each test to help the user determine
#' how many clusters is appropriate for their analysis.
#'
#' @param input_df A dataframe with variables as columns, samples as rows
#'
#' @return Prints 3 plots to test optimal cluster amounts for k-means clustering
#' @export
plot_optimal_cluster_tests <- function(input_df) {
    print(factoextra::fviz_nbclust(input_df, kmeans, "wss"))
    print(factoextra::fviz_nbclust(input_df, kmeans, "gap_stat"))
    print(factoextra::fviz_nbclust(input_df, kmeans, "silhouette"))
}

# Extracts clusters from kmeans_result object

#' Extract clusters from a kmeans_result object
#'
#' This function takes a kmeans_result object as input, and will extract
#' the cluster classifcation for each sample as a dataframe.
#' The assumption is that sample or patient IDs are stored as rownames
#' prior to clustering, for example at dimensionality reduction step.
#'
#' @param kmeans_cluster_result A result from the stats::kmeans() function
#'
#' @return A dataframe with columns "patient_id" and "cluster", where cluster
#' specifies clustering found via kmeans clustering.
#' @export
extract_clusters_from_kmeans <- function(kmeans_cluster_result) {
    # Extract clusters from kmeans cluster result
    clusters <- base::as.data.frame(kmeans_cluster_result[["cluster"]]) %>%
        tibble::rownames_to_column(var = "patient_id") %>%
        dplyr::rename("cluster" = dplyr::contains("cluster"))

    # Convert clusters to factor
    clusters[["cluster"]] <- base::as.factor(clusters[["cluster"]])

    return(clusters)
}

#' Perform kmeans clustering and extract cluster information
#'
#' This function performs kmeans clustering on a dataframe, where variables
#' are columns and samples are rows. It takes an argument, centers, which
#' should be the ideal amount of clustering, found through the
#' plot_optimal_cluster_tests function.
#'
#' @param random_seed An unsigned number that specifies a random seed.
#' This should bethe same random seed used in calculating the kmeans
#' clustering used in determining the optimal amount of clusters.
#' @param input_df A dataframe with variables as columns, samples as rows.
#' It is assumed the the rows are patient IDs.
#' @param centers The number of clusters to use in kmeans clustering.
#' @param nstart An unsigned number, defaults to 25. Determines how many
#' random sets should be chosen. ?kmeans for more information.
#'
#' @return A dataframe with columns of patient IDs and cluster information
perform_kmeans_and_extract_clusters <- function(
    random_seed,
    input_df,
    centers,
    nstart = 25
) {
    # Seed random seed
    set.seed(random_seed)

    # Perform kmeans
    kmeans_result <- stats::kmeans(input_df, centers = centers, nstart = nstart)

    # Printing cluster results
    print(factoextra::fviz_cluster(kmeans_result, data = input_df))

    result <- extract_clusters_from_kmeans(kmeans_result)

    return(result)
}