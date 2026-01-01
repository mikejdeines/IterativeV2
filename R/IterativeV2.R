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
  require(scCustomize)
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
  Idents(seurat.object) <- seurat.object$seurat_clusters
  seurat.object <- Rename_Clusters(seurat.object, seq(1, length(unique(seurat.object$seurat_clusters))))
  seurat.object$seurat_clusters <- Idents(seurat.object)
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
        cluster.object <- FindClusters(cluster.object, algorithm = 4, random.seed = 1)
        
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
CalculateDEScore <- function(seurat.object, cluster1, cluster2, pct.1, min.log2.fc){
  #' Calculates the DE score for a pair of clusters
  #' @param seurat.object a normalized, integrated Seurat object
  #' @param cluster1 the first cluster in the pair of clusters
  #' @param cluster2 the second cluster in the pair of clusters
  #' @param pct.1 the minimum gene expression fraction required for genes used to calculate the DE score
  #' @param min.log2.fc the minimum log2 fold change required for genes used to calculate the DE score
  #' @returns the DE score between the pair of clusters
  require(Seurat)
  if (sum(seurat.object$seurat_clusters == cluster1) < 3 || sum(seurat.object$seurat_clusters == cluster2) < 3) {
        return(0)
  }
  markers <- DGE.2samples(seurat.object, ident.1 = cluster1, ident.2 = cluster2, 
                          fc.thr = 1, min.pct = 0, max.pval = 1, min.count = 10,
                          icc = "i", df.correction = FALSE)
  markers <- markers[abs(markers$log2FC) < Inf, ]
  markers$p_val_adj <- p.adjust(markers$p.value, method = "BH")
  markers$pct.1 <- length(which(seurat.object@assays$RNA@counts[rownames(markers), colnames(seurat.object)[seurat.object$seurat_clusters == cluster1]] > 0)) / sum(seurat.object$seurat_clusters == cluster1)
  markers$pct.2 <- length(which(seurat.object@assays$RNA@counts[rownames(markers), colnames(seurat.object)[seurat.object$seurat_clusters == cluster2]] > 0)) / sum(seurat.object$seurat_clusters == cluster2)
  markers <- markers[abs(markers$log2FC) > min.log2.fc, ]
  markers <- markers[markers$pct.1 > pct.1 | markers$pct.2 > pct.1, ]
  markers$p_val_adj[markers$p_val_adj == 0] <- 1e-310
  markers$de_score <- -log10(markers$p_val_adj)
  markers$de_score[markers$de_score > 20] <- 20
  de_score <- sum(markers$de_score)
  return(de_score)
}
DGE.2samples <- function(object,features=NULL,ident.1=NULL,ident.2=NULL,
                         fc.thr=1,min.pct=0,max.pval=1,min.count=30,
                         icc="i",df.correction=FALSE) {
  IWT<-IterWghtTtest(object,features,ident.1,ident.2,fc.thr,min.pct,max.pval=max.pval,min.count,icc=icc,df.correction)
  Chi2<-Chi2Test(object,features=rownames(IWT),ident.1,ident.2,fc.thr=1,min.pct,max.pval=1,min.count)
  features<-intersect(rownames(Chi2),rownames(IWT))
  output<-IWT[features,]
  output[features,"Chi2.p.value"]<-Chi2[features,2]
  return(output)
}
IterWghtTtest <- function(object,features=NULL,ident.1=NULL,ident.2=NULL,fc.thr=1,min.pct=0,max.pval=1,min.count=30,
                          icc="i",df.correction=FALSE) {
  # Iterative weighted t-test for Seurat objects, see https://www.biorxiv.org/content/10.1101/2025.10.20.683496v1
  if (is.null(features)){
    GeneList <- rownames(object)
  } else {
    GeneList <- as.vector(features)
  }
  if (is.null(ident.1)|is.null(ident.2)) {stop("Two identities in a Seurat object must be defined using ident.1= and ident.2=")}
  output<-data.frame(log2FC=numeric(),p.value=numeric(), col3=numeric(),col4=numeric())                   # output dataframe
  colnames(output)<-c("log2FC","p.value", "Counts/Cell.1","Counts/Cell.2")                               # output dataframe
  object.1 <- subset(object, idents = ident.1)                                                           # ident.1 object
  Ci.1 <- as.matrix(object.1[["RNA"]]$counts)                                                            # ident.1 count matrix
  Ni.1 <- colSums(Ci.1)                                                                                  # ident.1 Ni vector
  Nc.1 <- ncol(Ci.1)                                                                                     # ident.1 number of cells
  Xi.1 <- Ci.1                                                                                           # normalized counts initiation     
  for (i in c(1:nrow(Ci.1))) {Xi.1[i,]=Ci.1[i,]/Ni.1}                                                    # ident.1 normalized counts
  object.2 <- subset(object, idents = ident.2)
  Ci.2 <- as.matrix(object.2[["RNA"]]$counts)
  Ni.2 <- colSums(Ci.2)
  Nc.2 <- ncol(Ci.2)
  Xi.2 <- Ci.2
  for (i in c(1:nrow(Ci.2))) {Xi.2[i,]=Ci.2[i,]/Ni.2}                                      
  rowi<-0
  print("Performing weighted t-test:")
  pb <- txtProgressBar(min = 0, max = length(GeneList), initial = 0, style = 3)                          # initialize progress bar# output row counter
  for (rownum in c(1:length(GeneList))) {                                                                # Main test loop
    AC.1 <- sum(Ci.1[GeneList[rownum],])                                                                 # ident.1 aggregated counts
    AC.2 <- sum(Ci.2[GeneList[rownum],])                                                                 # ident.2 aggregated counts
    Xi1<-Xi.1[GeneList[rownum],]                                                                         # ident.1 normalized counts
    Xi2<-Xi.2[GeneList[rownum],]                                                                         # ident.2 normalized counts
    if ((AC.1 >= min.count | AC.2 >= min.count)&                                                         # checking for minumum aggregated counts
        ((sum(Xi1!=0)/Nc.1>min.pct)|(sum(Xi2!=0)/Nc.2>min.pct))){                                        # checking for minumum expression
      wi.1<-ICCWeight(h=Ci.1[GeneList[rownum],],n=Ni.1,icc=icc)                                          # ident.1 weights
      wi.2<-ICCWeight(h=Ci.2[GeneList[rownum],],n=Ni.2,icc=icc)                                          # ident.2 weights
      fc <- sum(Xi1*wi.1)/sum(Xi2*wi.2)                                                                  # fold change
      if (!is.na(fc)){                                                                                   # removing 0/0 division
        if ((fc >= fc.thr | fc <= 1/fc.thr)&                                                             # checking FC threshold
            (sum(Xi1!=0)>=3|sum(Xi2!=0)>=3)) {                                                           # checking for at least 3 nonzero values in ident.1 or ident.2                                                      
          if(df.correction==TRUE) {wTtest <- as.numeric(alt.wttest2(Xi1,Xi2,wi.1,wi.2))}                 # weighted t-test with effective df
          else {wTtest <- as.numeric(alt.wttest(Xi1,Xi2,wi.1,wi.2))}                                     # weighted t-test without df correction
          if (wTtest <= max.pval) {                                                                      # checking p-value threshold    
            rowi<-rowi+1                                                                                 # output row counter
            output[[rowi,1]] <- as.numeric(format(log(fc,2), digits=4, scientific=FALSE))                # output
            output[[rowi,2]] <- as.numeric(format(wTtest, digits=4, scientific=FALSE))
            output[[rowi,3]] <- as.numeric(format(AC.1/Nc.1, digits=4, scientific=FALSE))
            output[[rowi,4]] <- as.numeric(format(AC.2/Nc.2, digits=4, scientific=FALSE))
            row.names(output)[rowi]<-GeneList[rownum]
          }  
        } 
      }
    }
    setTxtProgressBar(pb, rownum)                                                                        # progress bar update
  }
  close(pb)                                                                                              # close progress bar
  return(output)
}
Chi2Test <- function(object,features=NULL,ident.1=NULL,ident.2=NULL,fc.thr=1,min.pct=0,max.pval=1,min.count=30) {
  if (is.null(features)){
    GL <- rownames(object)
  } else {
    GL <- as.vector(features)
  }
  if (is.null(ident.1)|is.null(ident.2)) {stop("Two identities in a Seurat object must be defined using ident.1= and ident.2=")}
  output<-data.frame(log2FC=numeric(),p.value=numeric(),col3=numeric(),col4=numeric())                               # output dataframe
  colnames(output)<-c("log2FC","p.value","Counts/Cell.1","Counts/Cell.2")
  object.1 <- subset(object, idents = ident.1)                                                                       # object subsetting
  object.2 <- subset(object, idents = ident.2)
  Ci.1 <- as.matrix(object.1[["RNA"]]$counts)                                                                        # count matrices
  Ci.2 <- as.matrix(object.2[["RNA"]]$counts)
  Nc.1 <- ncol(Ci.1)                                                                                                 # number of cells
  Nc.2 <- ncol(Ci.2)
  AC.1 <- rowSums(Ci.1)                                                                                              # aggregate counts /gene
  AC.2 <- rowSums(Ci.2)
  TC.1 <- sum(AC.1)                                                                                                  # total counts /sample
  TC.2 <- sum(AC.2)
  rowi<-0
  print("Performing chi^2 test:")
  pb <- txtProgressBar(min = 0, max = length(GL), initial = 0, style = 3)                                            # initialize progress bar
  for (rownum in c(1:length(GL))) {                                                                                  # main loop, gene by gene
    if ((AC.1[GL[rownum]] >= min.count | AC.2[GL[rownum]] >= min.count)&                                             # checking for min.count and min.pct
        ((sum(Ci.1[GL[rownum],]!=0)/Nc.1>min.pct)|(sum(Ci.2[GL[rownum],]!=0)/Nc.2>min.pct))){
      ContTable <- matrix(c(TC.1 - AC.1[GL[rownum]], TC.2 - AC.2[GL[rownum]],                                        # contingency table for chi2 test
                            AC.1[GL[rownum]], AC.2[GL[rownum]]), nrow = 2, ncol = 2)
      fc <- as.numeric((AC.1[GL[rownum]]/TC.1)/(AC.2[GL[rownum]]/TC.2))                                              # fold change (FC)
      if(!is.na(fc)){                                                                                                # removing features with undefined FC (=0/0)
        if (fc >= fc.thr | fc <= 1/fc.thr) {                                                                         # checking FC threshold
          chisquare <- as.numeric(chisq.test(ContTable)$p.value)                                                     # chi2 test
          if (chisquare <= max.pval) {                                                                               # checking max.pval threshold
            rowi<-rowi+1
            output[[rowi,1]] <- as.numeric(format(log(fc,2), digits=4, scientific=FALSE))                            # output assembly
            output[[rowi,2]] <- as.numeric(format(chisquare, digits=4, scientific=FALSE))
            output[[rowi,3]] <- as.numeric(format(AC.1[GL[rownum]]/Nc.1, digits=4, scientific=FALSE))
            output[[rowi,4]] <- as.numeric(format(AC.2[GL[rownum]]/Nc.2, digits=4, scientific=FALSE))
            row.names(output)[rowi]<-GL[rownum]
          }  
        }
      }
    }
    setTxtProgressBar(pb, rownum)                                                                                    # progress bar update
  }
  close(pb)                                                                                                          # close progress bar
  return(output)
}
ICC.AN<-function(h,n){
  N=sum(n)                                                    #ANOVA ICC calculation
  k=length(n)
  n0=(1/(k-1))*(N-sum(n^2/N))
  MSw=(1/(N-k))*(sum(h)-sum(h^2/n))
  MSb=(1/(k-1))*(sum(h^2/n)-(1/N)*(sum(h))^2)
  if ((MSb+(n0-1)*MSw)==0) {ICC=0}                            #setting ICC=0 when ANOVA denominator = 0
  else ICC=(MSb-MSw)/(MSb+(n0-1)*MSw)
  if(ICC < 0) ICC=0                                           #resetting negative ICC to 0
  if(ICC > 1) ICC=1                                           #resetting ICC>1 to ICC=1
  return(ICC)
}
ICC.iter<-function(h,n){                                      #Iterative ICC calculation
  x=h/n
  w0<-n/sum(n)                                              # initial weights 
  x0<-sum(x*w0)                                             # initial weighted average count
  VarT0<-x0*(1-x0)/sum(n)                             # initial variance @ icc=0
  VarE0<-sum(w0^2*(x-x0)^2)/(1-sum(w0^2))                   # initial measured variance with w0 weights
  if (VarE0<=VarT0){icc=0}    
  else{
    f <- function(icc,h=h,n=n) {
      x = h/n                                                   #normalized counts
      wprop = n/(1 + icc*(n-1))                                 #proportional weights
      w <- wprop/sum(wprop)                                     #normalized weights
      x1 = sum(x*w)                                             #weighted average
      VarT<-x1*(1-x1)/(sum(wprop))                              #VarT 
      VarE<-sum(w^2*(x-x1)^2)/(1-sum(w^2))                      #VarE
      VarE-VarT                                                 #VarE-VarT
    }
    ur = uniroot(f, 0:1, n=n, h=h, check.conv=T, extendInt="downX", tol = 1e-4/max(n))
    icc = ur$root
    if(icc > 1) icc=1                                           #resetting ICC>1 to ICC=1
  }
  return(icc)
}
ICCWeight <- function(h,n,icc="i") {                     #Calculation of weights based on ICC
  Nc<-length(n)
  if(length(h)!=Nc) {stop("Unequal lengths of Ni and Nzi vectors")} 
  if(Nc<3) {stop("At least 3 cells are required in each ident")}
  if (sum(h!=0)<3){w<-rep(1/Nc,Nc)} #exit and return equal weights if the number cells/samples with nonzero counts is less than 3
  else {
    if(icc=="i"){icc=ICC.iter(h,n)} else {
      if(icc=="A"){icc=ICC.AN(h,n)} else {
        if (icc==0) {icc=0} else {
          if (icc==1) {icc=1} else {
            stop("Invalid icc, must be icc = 'i', 'A', 0, or 1")
          }
        }  
      }
    } 
    wprop = n/(1 + icc*(n-1))                                 #proportional weights
    w <- wprop/sum(wprop)
  }
  return(w)
}