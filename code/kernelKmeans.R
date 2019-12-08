## Kernel Kmeans function
# Perform Kernel Kmeans given a kernel matrix and a cluster number.
# The algorihtm randomly assigns a point to each cluster (center) and performs
# the reassignment iterations until convergence or max iterations.
# Optionally, the number of random center assigments can be specified so 
# that the algorithms repeat the previous process to each set an takes the best.
kernel.kmeans <- function(
  K, # Kernel matrix
  clusters, # Number of clusters
  nrandomSet = 5, # Number of random center start assignment to test
  iter.max=100, # Iteration maximum
  verbose=FALSE # Toggle verbose mode
) {
  N <- nrow(K)
  cnt <- lapply(1:nrandomSet, function(l) {sample(N, size=clusters, replace=FALSE)})
  withinss <- Inf
  for(i in 1:nrandomSet) {
    res <- .one.kkmeans(K,clusters,cnt[[i]],iter.max,verbose)
    new.withinss <- sum(res$withinss)
    if(withinss > new.withinss) {
      Z <- res$Z
      wss <- res$withinss
      withinss <- new.withinss
      centers <- res$centers
    }
  }
  return(Z)
}

.one.kkmeans <- function(
  K, # Kernel matrix
  clusters=clusters, # Number of clusters to consider
  centers=sample(1:nrow(K), size = clusters, replace = F), # First assignments
  iter.max=100, # Iteration maximum
  verbose=FALSE # Toggle verbose mode
) {
  # Assignment matrix
  Z <- matrix(0, nrow = nrow(K), ncol = clusters)
  # Assign centers to clusters
  for(i in 1:clusters) {
    Z[centers[i],i] <- 1
  }
  # First iteration with given centers
  res <- .assign(K,Z,clusters)
  Z <- res$Z
  # Start iterations
  iter <- 1
  num.assign <- nrow(K)
  # Stop if num.assign are 0 (or if iter.max reached)
  while(iter < iter.max & num.assign != 0) {
    iter <- iter+1
    res <- .assign(K,Z,clusters)
    Z <- res$Z
    num.assign <- res$num.assign        
  }
  # Compute centers and witin cluster sum of squares
  clust <- apply(Z,1,which.max)
  cent <- lapply(1:clusters, function(c) .getCenters(c, clust, K))
  withinss <- unlist(lapply(1:clusters, function(c) .getWithinSS(c,clust,cent,K)))
  if(verbose) {
    print(paste("Finished Kernel Kmeans iterations with",iter,
      "iterations."))
  }
  return(list(Z=Z, withinss=withinss, centers=cent))
}

.assign <- function(
  K,
  Z,
  clusters
) {
  # Cluster size matrix
  L <- diag(apply(Z, 2, function(x) 1/sum(x)))
  # Within cluster kernel sum
  g <- rep(0, clusters)
  for(l in 1:clusters) {
    G <- K[which(Z[,l] > 0),which(Z[,l] > 0)]
    g[l] <- sum(G)*(L[l,l]^2)
  }
  # Compute point-cluster kernels
  F <- K %*% Z %*% (-2*L)
  # Sum within cluster and point-cluster kernels
  d <- apply(F, 1, function(x) x + g)
  # Get new assignments
  num.assign <- 0
  rassingn <- apply(d, 2, which.min)
  for(i in 1:nrow(K)) {
    oldassign <- ifelse(length(which(Z[i,]>0))==0,Inf,which(Z[i,]>0))
    if(oldassign != rassingn[i]) {
      num.assign <- num.assign+1
      if(oldassign!=Inf) {
        Z[i,oldassign] <- 0
      }
      Z[i,rassingn[i]] <- 1
    }
  }
  return(list(Z=Z,
              num.assign=num.assign))
}

.getCenters <- function(
  c,
  clust,
  K
) {
  G <- K[which(clust==c),which(clust==c),drop=FALSE]
  if(nrow(G)>1) {
    cent <- colMeans(G)
  } else {
    cent <- G[,,drop=TRUE]
  }
  return(cent)
}

.getWithinSS <- function(
  c,
  clust,
  cent,
  K
) {
  G <- K[which(clust==c),which(clust==c)]
  SS <- G - cent[[c]]
  wss <- sum(SS^2)
  return(wss)
}
