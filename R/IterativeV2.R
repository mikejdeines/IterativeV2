IterativeClustering <- function(seurat.object, max.iterations = 20,
                                min.cluster.size = 10, min.de.score = 150, pct.1 = 0.5,
                                min.log2.fc = 2, n.dims = 30, integration_method = "None", batch_key = "None", n.cores = 1) {
  #' Performs iterative Leiden clustering on a Seurat object
  #' @param seurat.object a normalized, integrated Seurat object
  #' @param max.iterations the maximum number of iterations to perform. Default 20.
  #' @param min.cluster.size the minimum number of cells in a cluster. Default 10.
  #' @param min.de.score the minimum DE score required to consider clusters biologically-distinct. Default 150.
  #' @param pct.1 the minimum gene expression fraction required for genes used to calculate the DE score. Default 0.5.
  #' @param min.log2.fc the minimum log2 fold change required for genes used to calculate the DE score. Default 2.
  #' @param n.dims number of dimensions used to create the sNN graph and calculate centroids. Default 30.
  #' @param integration_method the integration method used ("None", "Harmony", "CCA", or "RPCA"). Default "None".
  #' @param batch_key the metadata key for batch information. Default "None".
  #' @param n.cores number of cores to use for DGE.2samples. Default 1.
  #' @returns a Seurat object with iterative clustering results in the seurat_clusters metadata value
  require(Seurat)
  require(scCustomize)
  iteration <- 1
  prev.n.clusters <- 1
  seurat.object$seurat_clusters <- 1
  repeat {
    seurat.object <- RunClusteringIteration(seurat.object, min.cluster.size, min.de.score, pct.1, min.log2.fc, n.dims, integration_method, batch_key, n.cores)
    curr.n.clusters <- length(unique(seurat.object$seurat_clusters))
    message(paste0("Iteration ", iteration, ": ", curr.n.clusters, " clusters identified."))
    if (iteration >= max.iterations || curr.n.clusters == prev.n.clusters) {
      break
    }
    prev.n.clusters <- curr.n.clusters
    iteration <- iteration + 1
  }
  Idents(seurat.object) <- seurat.object$seurat_clusters
  seurat.object <- Rename_Clusters(seurat.object, seq(1, length(unique(seurat.object$seurat_clusters))))
  seurat.object$renamed_clusters <- Idents(seurat.object)
  return(seurat.object)
}
RunClusteringIteration <- function(seurat.object, min.cluster.size, min.de.score, pct.1, min.log2.fc, n.dims, integration_method, batch_key, n.cores){
  #' Runs a single Leiden clustering iteration
  #' @param seurat.object a normalized, integrated Seurat object
  #' @param min.cluster.size the minimum number of cells in a cluster
  #' @param min.de.score the minimum DE score required to consider clusters biologically-distinct
  #' @param pct.1 the minimum gene expression fraction required for genes used to calculate the DE score
  #' @param min.log2.fc the minimum log2 fold change required for genes used to calculate the DE score
  #' @param n.dims number of dimensions used to create the sNN graph and calculate centroids
  #' @param integration_method the integration method used ("None", "Harmony", "CCA", or "RPCA")
  #' @param batch_key the metadata key for batch information
  #' @param n.cores number of cores to use for DGE.2samples
  #' @returns a Seurat object with the results of the iteration in the seurat_clusters metadata value
  require(Seurat)
    unique_clusters <- unique(seurat.object$seurat_clusters)
    
    for (cluster_id in unique_clusters){
        Idents(seurat.object) <- seurat.object$seurat_clusters
        cluster.object <- subset(seurat.object, idents = cluster_id)
        if (ncol(cluster.object) < min.cluster.size) {
            next
        }
        if (integration_method == "None" || nlevels(factor(cluster.object[[batch_key]][,1])) < 2){
          dim.reduction <- "pca"
          message("Finding variable features and scaling data...")
          cluster.object <- FindVariableFeatures(cluster.object, verbose = FALSE)
          cluster.object <- ScaleData(cluster.object, verbose = FALSE, features = VariableFeatures(cluster.object))
          message("Running PCA...")
          cluster.object <- RunPCA(cluster.object, npcs = n.dims, verbose = FALSE, features = VariableFeatures(cluster.object), approx = TRUE)
        } else if (integration_method == "Harmony" && nlevels(factor(cluster.object[[batch_key]][,1])) >= 2){
          dim.reduction <- "harmony"
          message("Finding variable features and scaling data...")
          cluster.object <- FindVariableFeatures(cluster.object, verbose = FALSE)
          cluster.object <- ScaleData(cluster.object, verbose = FALSE, features = VariableFeatures(cluster.object))
          message("Running PCA...")
          cluster.object <- RunPCA(cluster.object, npcs = n.dims, verbose = FALSE, features = VariableFeatures(cluster.object), approx = TRUE)
          cluster.object[["RNA"]] <- split(cluster.object[["RNA"]], f = cluster.object[[batch_key]][,1])
          message("Running Harmony integration...")
          cluster.object <- IntegrateLayers(cluster.object, method = HarmonyIntegration)
          message("Joining integrated layers...")
          cluster.object[["RNA"]] <- JoinLayers(cluster.object[["RNA"]])
        } else if (integration_method == "CCA" && nlevels(factor(cluster.object[[batch_key]][,1])) >= 2){
          dim.reduction <- "integrated.dr"
          message("Finding variable features and scaling data...")
          cluster.object <- FindVariableFeatures(cluster.object, verbose = FALSE)
          cluster.object <- ScaleData(cluster.object, verbose = FALSE, features = VariableFeatures(cluster.object))
          message("Running PCA...")
          cluster.object <- RunPCA(cluster.object, npcs = n.dims, verbose = FALSE, features = VariableFeatures(cluster.object), approx = TRUE)
          cluster.object[["RNA"]] <- split(cluster.object[["RNA"]], f = cluster.object[[batch_key]][,1])
          message("Running CCA integration...")
          cluster.object <- IntegrateLayers(cluster.object, method = CCAIntegration)
          message("Joining integrated layers...")
          cluster.object[["RNA"]] <- JoinLayers(cluster.object[["RNA"]])
        } else if (integration_method == "RPCA" && nlevels(factor(cluster.object[[batch_key]][,1])) >= 2){
          dim.reduction <- "integrated.dr"
          message("Finding variable features and scaling data...")
          cluster.object <- FindVariableFeatures(cluster.object, verbose = FALSE)
          cluster.object <- ScaleData(cluster.object, verbose = FALSE, features = VariableFeatures(cluster.object))
          message("Running PCA...")
          cluster.object <- RunPCA(cluster.object, npcs = n.dims, verbose = FALSE, features = VariableFeatures(cluster.object), approx = TRUE)
          cluster.object[["RNA"]] <- split(cluster.object[["RNA"]], f = cluster.object[[batch_key]][,1])
          message("Running RPCA integration...")
          cluster.object <- IntegrateLayers(cluster.object, method = RPCAIntegration)
          message("Joining integrated layers...")
          cluster.object[["RNA"]] <- JoinLayers(cluster.object[["RNA"]])
        } else {
          stop("Invalid integration method or batch key.")
        }
        message("Creating sNN graph...")
        if (ncol(cluster.object) < 20){
          cluster.object <- FindNeighbors(cluster.object, dims = 1:n.dims, reduction = dim.reduction, k.param = floor(ncol(cluster.object)/2), verbose = FALSE)
        }
        else {
          cluster.object <- FindNeighbors(cluster.object, dims = 1:n.dims, reduction = dim.reduction, verbose = FALSE)
        }
        message("Finding clusters with Leiden algorithm...")
        cluster.object <- FindClusters(cluster.object, algorithm = 4, random.seed = 1, leiden_method = "igraph", resolution = 1)
        message(paste(length(unique(cluster.object$seurat_clusters)), "sub-clusters found. Merging similar clusters..."))
        # Initialize for merge loop
        sub_clusters <- unique(cluster.object$seurat_clusters)
        if (length(sub_clusters) > 1) {
            centroids <- FindCentroids(cluster.object, n.dims, dim.reduction)
            dist_matrix <- as.matrix(dist(centroids))
        }
        
        repeat {
                sub_clusters <- unique(cluster.object$seurat_clusters)
                if (length(sub_clusters) <= 1) {
                    break
                }
                merged_pairs <- list()
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
                    de_score <- CalculateDEScore(cluster.object, cluster1, cluster2, pct.1, min.log2.fc, n.cores)
                    if (de_score < min.de.score || sum(cluster.object$seurat_clusters == cluster2) < min.cluster.size || sum(cluster.object$seurat_clusters == cluster1) < min.cluster.size) {
                        cluster.object$seurat_clusters[cluster.object$seurat_clusters == cluster2] <- cluster1
                        merged_pairs <- c(merged_pairs, pair_key)
                        merged <- TRUE
                        # Recalculate centroids and distance matrix after merge
                        centroids <- FindCentroids(cluster.object, n.dims, dim.reduction)
                        dist_matrix <- as.matrix(dist(centroids))
                        break
                    }
                }
                
                # If no merges happened this iteration, exit
                if (!merged) {
                    break
                }
        }
        
        cluster.cells <- WhichCells(cluster.object)
        sub_clusters_final <- unique(cluster.object$seurat_clusters)
        if (length(sub_clusters_final) == 1) {
            # If only one subcluster, label as parent_1
            seurat.object$seurat_clusters[cluster.cells] <- paste0(cluster_id, "_1")
        } else {
            for (i in seq_along(sub_clusters_final)) {
                k <- sub_clusters_final[i]
                sub_cluster_cells <- WhichCells(cluster.object, idents = k)
                seurat.object$seurat_clusters[sub_cluster_cells] <- paste0(cluster_id, "_", i)
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
CalculateDEScore <- function(seurat.object, cluster1, cluster2, pct.1, min.log2.fc, n.cores){
  #' Calculates the DE score for a pair of clusters
  #' @param seurat.object a normalized, integrated Seurat object
  #' @param cluster1 the first cluster in the pair of clusters
  #' @param cluster2 the second cluster in the pair of clusters
  #' @param pct.1 the minimum gene expression fraction required for genes used to calculate the DE score
  #' @param min.log2.fc the minimum log2 fold change required for genes used to calculate the DE score
  #' @param n.cores number of cores to use for DGE.2samples
  #' @returns the DE score between the pair of clusters
  require(Seurat)
  if (sum(seurat.object$seurat_clusters == cluster1) < 3 || sum(seurat.object$seurat_clusters == cluster2) < 3) {
        return(0)
  }
  markers <- DGE.2samples(seurat.object, ident.1 = cluster1, ident.2 = cluster2, 
                          fc.thr = 1, min.pct = 0, max.pval = 1, min.count = 10,
                          icc = "i", df.correction = FALSE, n.cores = n.cores)
  markers <- markers[abs(markers$log2FC) < Inf, ]
  
  if (nrow(markers) == 0) {
    return(0)
  }
  
  markers$p_val_adj <- p.adjust(markers$p.value, method = "BH")
  
  # Get count matrix
  counts_matrix <- seurat.object@assays[["RNA"]]@layers$counts
  cells_cluster1 <- which(seurat.object$seurat_clusters == cluster1)
  cells_cluster2 <- which(seurat.object$seurat_clusters == cluster2)
  
  # Calculate pct.1 and pct.2 for each gene using vectorized operations
  gene_indices <- match(rownames(markers), rownames(seurat.object))
  counts_cluster1 <- counts_matrix[gene_indices, cells_cluster1, drop = FALSE]
  counts_cluster2 <- counts_matrix[gene_indices, cells_cluster2, drop = FALSE]
  
  # Use Matrix::rowSums for sparse matrices
  if (inherits(counts_matrix, "dgCMatrix") || inherits(counts_matrix, "IterableMatrix")) {
    markers$pct.1 <- Matrix::rowSums(counts_cluster1 > 0) / length(cells_cluster1)
    markers$pct.2 <- Matrix::rowSums(counts_cluster2 > 0) / length(cells_cluster2)
  } else {
    markers$pct.1 <- rowSums(counts_cluster1 > 0) / length(cells_cluster1)
    markers$pct.2 <- rowSums(counts_cluster2 > 0) / length(cells_cluster2)
  }
  
  markers <- markers[abs(markers$log2FC) > min.log2.fc, ]
  markers <- markers[markers$pct.1 > pct.1 | markers$pct.2 > pct.1, ]
  markers$p_val_adj[markers$p_val_adj == 0] <- 1e-310
  markers$de_score <- -log10(markers$p_val_adj)
  markers$de_score[markers$de_score > 20] <- 20
  de_score <- sum(markers$de_score)
  return(de_score)
}