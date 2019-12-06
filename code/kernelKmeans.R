## Kernel Kmeans function
# Perform Kernel Kmeans given a kernel matrix and a cluster number.
# Optionally, the initial assignments can be specified. If not, the
# function will randomly assign instances to all classes.
kernel.kmeans <- function(
  K, # Kernel matrix
  clusters, # Number of clusters
  assignments=NA, # Assignment matrix nrow(K) x clusters with 0-1
  iter.max=100 # Iteration maximum
) {
  # Assignment matrix
  if((!is.null(dim(assignments))
     && dim(assignments)[1] == nrow(K)
     && dim(assignments)[2] == clusters)) {
    Z <- assignments
  } else {
    # Equivalent assignments
    Z <- matrix(0, nrow = nrow(K), ncol = clusters)
    part <- nrow(K)/clusters +1
    rassign <- sapply(1:nrow(K)/part, function(x) floor(x)+1)
    for(i in 1:nrow(K)) {
      Z[i,rassign[i]] <- 1
    }
  }  
  # Cluster size matrix
  L <- diag(apply(Z, 2, function(x) 1/sum(x)))
  # Within cluster kernel sum
  g <- rep(0, clusters)
  for(l in 1:clusters) {
    G <- K[which(Z[,l] > 0),which(Z[,l] > 0)]
    g[l] <- sum(G)/L[l,l]^-2
  }  
  # Start iteration
  iter <- 1
  num.assign <- nrow(K)
  # Stop if num.assign are 0 (or if iter.max reached)
  while(iter < iter.max & num.assign != 0) {
    iter <- iter+1
    num.assign <- 0
    # Compute point-cluster kernels
    F <- K %*% Z %*% (-2*L)
    # Sum within cluster and point-cluster kernels
    d <- apply(F, 1, function(x) x + g)
    # Get new assignments
    rassingn <- apply(d, 2, which.min)
    for(i in 1:nrow(K)) {
      oldassign <- which(Z[i,] > 0)
      if(oldassign != rassingn[i]) {
        num.assign <- num.assign+1
        Z[i,oldassign] <- 0
        Z[i,rassingn[i]] <- 1
      }
    }
    # Recompute cluster size matrix
    L <- diag(apply(Z, 2, function(x) 1/sum(x)))    
    # Recompute within cluster kernels
    for(l in 1:clusters) {
      G <- K[which(Z[,l] > 0),which(Z[,l] > 0)]
      g[l] <- sum(G)/L[l,l]^-2
    }
  }
  print(paste("Finished Kernel Kmeans iterations with",iter,
              "iterations."))
  return(Z)
}
