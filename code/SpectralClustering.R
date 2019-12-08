
spectralClustering <- function(
  S, #similarity matrix, S,
  clusters, # Number of clusters
  laplacian = "normalizedRandomWalk", # or unnormalized,normalizedSymmetric
  nrandomSet 
  ){
  # degree matrix D where each diagonal value is the degree of the respective vertex 
  # and all other positions are zero
  D <- diag(apply(S, 1, sum))
  # Calculate the Laplacian
  # which laplacian should we pick? look at the degree distribution of the similarity graph
  if(laplacian=="unnormalized"){ # Shi and Malik (2000)
    L <- D - S
  } else if(laplacian == "normalizedSymmetric"){ #Ng, Jordan, and Weiss (2002)
    L <- .sqrtm(D) %*% S %*% .sqrtm(D)
    # additional multiplication D might lead to undesired artifacts
  } else if(laplacian=="normalizedRandomWalk"){ # preferred by Von Luxburg, U. (2007). 
    # cluster indicator vectors 1_{A_i}
    L <- .ginv(D) %*% S
  }
  # Compute eigenvectors and eigenvalues of L
  evL <- eigen(L, symmetric=TRUE)
  # choose eigenvectors
  n <- ncol(S)
  Z <- evL$vectors[,((n-clusters+1):n)]
  # run K-means algorithm on that space to get k clusters 
  C <- kmeans(Z, centers=clusters,nstart = nrandomSet, iter.max =100) 
  classes <- C$cluster
  assignments <- matrix(0, nrow = length(classes), ncol = clusters)
  for(i in 1:length(classes)) {
    assignments[i,classes[i]] <- 1
  }
  return(assignments)
}


.ginv <- function (X, tol = sqrt(.Machine$double.eps))
{
  if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X))
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
  if (all(Positive))
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive))
    array(0, dim(X)[2:1])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}

.sqrtm <- function(x)
{
  tmpres <- eigen(x)
  V <- t(tmpres$vectors)
  D <- tmpres$values
  if(is.complex(D))
    D <- Re(D)
  D <- pmax(D,0)
  return(crossprod(V*sqrt(D),V))
}