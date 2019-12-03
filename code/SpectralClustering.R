
######################################################
##############      functions ##############
######################################################

# A Gaussian kernel, s, to calculate the similarity between two points.
s <- function(x1, x2, alpha=1) {
  exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

#compute a restricted (or filtered) “affinity” between vertices using k-nearest neighbors.
make.affinity <- function(S, n.neighboors=2) {
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
  A
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
######################################################
##############      algorithm ##############
######################################################
# 1.Compute the similarity matrix, S, between all points 
# 2.Calculate the affinity matrix, A, from S by applying k-nearest neighbors algorithm
# 3.Compute the weighted degree diagonal matrix, D, from A by summing across each row
# 4.Calculate the unnormalized Laplacian, U, by subtracting A from D
# 5.Compute eigenvectors and eigenvalues of U.
# 6.Perform k-means clustering on k smallest eigenvalues, ignoring the smallest (constant) eigenvector
######################################################
library(kernlab)
data(spirals)

# 1.Compute the similarity matrix, S, between all points 
S <- make.similarity(spirals, s)
# 2.Calculate the affinity matrix, A, from S by applying k-nearest neighbors algorithm
# n.neighboors parameter needs to be adjusted depending on the data 
A <- make.affinity(S, 3)  
# 3. Compute the weighted degree diagonal matrix, D, from A by summing across each row
D <- diag(apply(A, 1, sum))
# 4.Calculate the unnormalized Laplacian, U, by subtracting A from D
##  
# i.Simple Laplacian L=D-A
# ii.Normalized Laplacian L_{N}=D^{-1/2}LD^{-1/2}
# iii.Generalized Laplacian L_{G} = D^{-1}L
# iv.Relaxed Laplacian L_{\rho} = L-\rho D 
# v.Ng, Jordan, & Weiss Laplacian L_{NJW}=D^{-1/2}AD^{-1/2}, and where A_{i,i}=0 
# vi.smoothed Kernel for Kmeans Kernel Clustering  K=\sigma D^{-1}+D^{-1}AD^{-1} 
##
#Left multiplying by a diagonal matrix is akin to scaling the rows, 
#but right multiplying by a diagonal matrix is akin to scaling the columns.  
#The generalized Laplacian results in a right-stochastic Markov matrix ; the normalized Laplacian does not.
##  
U <- D - A
# 5.Compute eigenvectors and eigenvalues of U.
evL <- eigen(U, symmetric=TRUE)
# performing k-means on eigenvectors of a similarity kernel applied to the original data.
k   <- 2
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
km <- kmeans(Z, centers=k, nstart=5)

