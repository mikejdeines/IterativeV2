IterativeClustering <- function(seurat.object, max.iterations = 20,
                                min.cluster.size = 10, min.de.score = 150, pct.1 = 0.5,
                                min.log2.fc = 2, n.dims = 30, dim.reduction = "pca") {
  #' Performs iterative Leiden clustering on a Seurat object
  #' @param seurat.object a normalized, integrated Seurat object
  #' @param max.iterations the maximum number of iterations to perform. Default 20.
  #' @param min.cluster.size the minimum number of cells in a cluster. Default 10.
  #' @param min.de.score the minimum DE score required to consider clusters biologically-distinct. Default 150.
  #' @param pct.1 the minimum gene expression fraction required for genes used to calculate the DE score. Default 0.5.
  #' @param min.log2.fc the minimum log2 fold change required for genes used to calculate the DE score. Default 2.
  #' @param n.dims number of dimensions used to create the sNN graph and calculate centroids. Default 30.
  #' @param dim.reduction the dimensional reduction used to create the sNN graph and calculate centroids. Default pca.
  #' @returns a Seurat object with iterative clustering results in the seurat_clusters metadata value
  require(Seurat)
  iteration <- 1
  prev.n.clusters <- 1
  seurat.object$seurat_clusters <- 1
  repeat {
    seurat.object <- RunClusteringIteration(seurat.object, min.cluster.size, min.de.score, pct.1, min.log2.fc, n.dims, dim.reduction)
    curr.n.clusters <- length(unique(seurat.object$seurat_clusters))
    if (iteration >= max.iterations || curr.n.clusters == prev.n.clusters) {
      break
    }
    prev.n.clusters <- curr.n.clusters
    iteration <- iteration + 1
  }
  return(seurat.object)
}
RunClusteringIteration <- function(seurat.object, min.cluster.size, min.de.score, pct.1, min.log2.fc, n.dims, dim.reduction){
  #' Runs a single Leiden clustering iteration
  #' @param seurat.object a normalized, integrated Seurat object
  #' @param min.cluster.size the minimum number of cells in a cluster
  #' @param min.de.score the minimum DE score required to consider clusters biologically-distinct
  #' @param pct.1 the minimum gene expression fraction required for genes used to calculate the DE score
  #' @param min.log2.fc the minimum log2 fold change required for genes used to calculate the DE score
  #' @param n.dims number of dimensions used to create the sNN graph and calculate centroids
  #' @param dim.reduction the dimensional reduction used to create the sNN graph and calculate centroids
  #' @returns a Seurat object with the results of the iteration in the seurat_clusters metadata value
  require(Seurat)
  unique_clusters <- unique(seurat.object$seurat_clusters)
    
  for (cluster_id in unique_clusters){
        Idents(seurat.object) <- seurat.object$seurat_clusters
        cluster.object <- subset(seurat.object, idents = cluster_id)
        if (ncol(cluster.object) < min.cluster.size) {
            next
        }
        if (ncol(cluster.object) < 20){
          cluster.object <- FindNeighbors(cluster.object, dims = 1:n.dims, reduction = dim.reduction, k.param = floor(ncol(cluster.object)/2), verbose = FALSE)
        }
        else {
          cluster.object <- FindNeighbors(cluster.object, dims = 1:n.dims, reduction = dim.reduction, verbose = FALSE)
        }
        cluster.object <- FindClusters(cluster.object, resolution = 0.8, algorithm = 4, method = "igraph")
        
        repeat {
                sub_clusters <- unique(cluster.object$seurat_clusters)
                prev_n_subclusters <- length(sub_clusters)
                if (prev_n_subclusters <= 1) {
                    break
                }
                changed <- FALSE
                merged_pairs <- list()
                centroids <- FindCentroids(cluster.object, n.dims, dim.reduction)
                dist_matrix <- as.matrix(dist(centroids))
                cluster_map <- setNames(seq_along(sub_clusters), sub_clusters)
                merged <- FALSE
                for (cluster1 in sub_clusters) {
                    if (!(as.character(cluster1) %in% names(cluster_map))) {
                        next
                    }
                    j <- cluster_map[[as.character(cluster1)]]
                    distances <- dist_matrix[j, ]
                    distances[j] <- Inf
                    closest_idx <- which.min(distances)
                    cluster2 <- sub_clusters[closest_idx]
                    pair_key <- paste(sort(c(cluster1, cluster2)), collapse = "_")
                    if (pair_key %in% merged_pairs) {
                        next
                    }
                    de_score <- CalculateDEScore(cluster.object, cluster1, cluster2, pct.1, min.log2.fc)
                    if (de_score < min.de.score || sum(cluster.object$seurat_clusters == cluster2) < min.cluster.size || sum(cluster.object$seurat_clusters == cluster1) < min.cluster.size) {
                        cluster.object$seurat_clusters[cluster.object$seurat_clusters == cluster2] <- cluster1
                        changed <- TRUE
                        merged_pairs <- c(merged_pairs, pair_key)
                        merged <- TRUE
                        break
                    }
                }
                if (merged) {
                    next
                }
                sub_clusters_new <- unique(cluster.object$seurat_clusters)
                if (length(sub_clusters_new) == prev_n_subclusters || !changed) {
                    break
                }
        }
        
        cluster.cells <- WhichCells(cluster.object)
        sub_clusters_final <- unique(cluster.object$seurat_clusters)
        if (length(sub_clusters_final) == 1) {
            seurat.object$seurat_clusters[cluster.cells] <- cluster_id
        } else {
            for (k in sub_clusters_final) {
                sub_cluster_cells <- WhichCells(cluster.object, idents = k)
                seurat.object$seurat_clusters[sub_cluster_cells] <- paste0(cluster_id, "_", k)
            }
        }
    }
    return(seurat.object)
  }
FindCentroids <- function(seurat.object, n.dims, dim.reduction) {
  #' Calculates centroids for each cluster in dim.reduction space
  #' @param seurat.object a normalized, integrated Seurat object
  #' @param n.dims number of dimensions used to calculate centroids
  #' @param dim.reduction the dimensional reduction used to calculate centroids
  #' @returns centroids for each cluster in dim.reduction space
  require(Seurat)
  embeddings <- Embeddings(seurat.object, reduction = dim.reduction)[, 1:n.dims]
  clusters <- seurat.object$seurat_clusters
  centroids <- aggregate(embeddings, by = list(cluster = clusters), FUN = mean)
  rownames(centroids) <- centroids$cluster
  centroids$cluster <- NULL
  return(centroids)
}
CalculateDEScore <- function(seurat.object, cluster1, cluster2, pct.1, min.log2.fc){
  #' Calculates the DE score for a pair of clusters
  #' @param seurat.object a normalized, integrated Seurat object
  #' @param cluster1 the first cluster in the pair of clusters
  #' @param cluster2 the second cluster in the pair of clusters
  #' @param pct.1 the minimum gene expression fraction required for genes used to calculate the DE score
  #' @param min.log2.fc the minimum log2 fold change required for genes used to calculate the DE score
  #' @returns the DE score between the pair of clusters
  require(Seurat)
  markers <- FindMarkers(seurat.object, ident.1 = cluster1, ident.2 = cluster2, min.pct = 0.1, logfc.threshold = min.log2.fc, recorrect_umi = FALSE)
  markers <- markers[abs(markers$avg_log2FC) > min.log2.fc, ]
  markers <- markers[markers$pct.1 > pct.1 | markers$pct.2 > pct.1, ]
  markers$p_val_adj[markers$p_val_adj == 0] <- 1e-310
  markers$de_score <- -log10(markers$p_val_adj)
  markers$de_score[markers$de_score > 20] <- 20
  de_score <- sum(markers$de_score)
  return(de_score)
}
