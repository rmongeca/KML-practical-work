
spectralClustering <- function(
  S, #similarity matrix, S,
  n.neighboors=round(log(dim(S)[1])), #choose k as per the assymptotic connectivity results
  clusters, # Number of clusters
  laplacian = "normalizedRandomWalk", # or unnormalized,normalizedSymmetric
  method = "kmeans" #pam or discretization
  ){
  # transform a given set of x_1,...,x_n with pairwise similarities S into a graph
  # connect vertex v_i with v_j if v_j is among the k-nearest neighboor of v_i
  # leads to a directed graph since neighboorhood relationship is not symmetric
  N <- length(S[,1])
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    # A must be made of positive values and be symmetric.
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        # solution applied to make the graph undirected so A is symmetric 
        # if A_ij is selected as a nearest neighboor, so will A_ji.
        A[j,i] <- S[i,j] 
      }
    }
  }
  # degree matrix D where each diagonal value is the degree of the respective vertex 
  # and all other positions are zero
  D <- diag(apply(A, 1, sum))
  # Calculate the Laplacian,
  # matrix power operator: computes M^power (M must be diagonalizable)
  "%^%" <- function(M, power)
    with(eigen(M), vectors %*% (values^power * solve(vectors)))
  # which laplacian ? look at the degree distribution of the similarity graph
  if(laplacian=="unnormalized"){ # Shi and Malik (2000)
    L <- D - A
  } else if(laplacian == "normalizedSymmetric"){ #Ng, Jordan, and Weiss (2002)
    L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))
    # additional multiplication D might lead to undesired
  } else if(laplacian=="normalizedRandomWalk"){ # preferred by Von Luxburg, U. (2007). 
    # cluster indicator vectors 1_{A_i}
    L <- (D %^% (-1)) %*% A
  }
  # Compute eigenvectors and eigenvalues of L
  # since we use the k nearest neighbor graph obtained Laplacian will be sparse
  evL <- eigen(L, symmetric=TRUE)
  # choose eigenvectors
  n <- ncol(S)
  Z <- evL$vectors[,((n-clusters+1):n)]
  if(method=="kmeans"){
    C<- kmeans(Z, centers=clusters, nstart=5)
    classes <- C$cluster
  }else if(method=="pam"){
    classes <- cluster::pam(x = Z, k = clusters,cluster.only=TRUE,diss = FALSE)
  }else if(method=="discretization"){
    res <- sort(abs(evL$values),index.return = TRUE)
    U <- evL$vectors[,res$ix[1:clusters]]
    normalize <- function(x) x / sqrt(sum(x^2))
    if(laplacian == "normalizedSymmetric"){
      U = t(apply(U,1,normalize))
    }
    eigDiscrete = .discretisation(U)
    eigDiscrete = eigDiscrete$discrete
    classes = apply(eigDiscrete,1,which.max)
  }
  assignments <- matrix(0, nrow = length(classes), ncol = length(unique(classes)))
  for(i in 1:length(classes)) {
    assignments[i,classes[i]] <- 1
  }
  return(assignments)
}


######## Internal Functions ##############

.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}

