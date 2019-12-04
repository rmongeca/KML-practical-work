
############## FUNCTION ###################

spectralClustering <- function(
  S, #similarity matrix, S,
  n.neighboors=2, #k-nearest neighbors algorithm
  clusters, # Number of clusters
  laplacian = "unnormalized" # or normalized
  ){
  #Calculate the affinity matrix, A, from S by applying k-nearest neighbors algorithm
  N <- length(S[,1])
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  # Compute the weighted degree diagonal matrix, D, from A by summing across each row
  D <- diag(apply(A, 1, sum))
  # Calculate the Laplacian,
  if(laplacian=="unnormalized"){
    L <- D - A
  } else if(laplacian == "normalized"){
    L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2)) 
  }
  # Compute eigenvectors and eigenvalues of L
  evL <- eigen(L, symmetric=TRUE)
  # find the k smallest eigenvectors (ignoring the trivial constant eigenvector)
  Z   <- evL$vectors[,(ncol(evL$vectors)-clusters+1):ncol(evL$vectors)]
  # performing k-means on eigenvectors of a similarity kernel 
  km <- kmeans(Z, centers=clusters, nstart=5)
  assignments <- km$cluster
  return(assignments)
}


############## TEST ###################

library(kernlab)
data(spirals)
data <- spirals
rm(spirals)

# A Gaussian kernel, s, to calculate the similarity between two points.
s <- function(x1, x2, alpha=1) {
  exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

make.similarity <- function(my.data, similarity) {
  N <- nrow(my.data)
  S <- matrix(rep(NA,N^2), ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      if (i!=j) {
        S[i,j] <- similarity(my.data[i,], my.data[j,])
      } else {
        S[i,j] <- 0
      }
    }
  }
  S
}
# similarity matrix, S,  
S <- make.similarity(data, s)

spectralCl <- spectralClustering(S = S,
                                 n.neighboors=3,
                                 clusters = 2, 
                                 laplacian = "unnormalized")
