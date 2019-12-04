library(kernlab)

#######  Mock data ###### 
make.sinusoidals <- function(m,noise=0.2) 
{
  x1 <- c(1:2*m)
  x2 <- c(1:2*m)
  for (i in 1:m) {
    x1[i] <- (i/m) * pi
    x2[i] <- sin(x1[i]) + rnorm(1,0,noise)
  }
  for (j in 1:m) {
    x1[m+j] <- (j/m + 1/2) * pi
    x2[m+j] <- cos(x1[m+j]) + rnorm(1,0,noise)
  }
  target <- as.factor(c(rep(+1,m),rep(-1,m)))
  return(data.frame(x1,x2,target))
}
data <- make.sinusoidals(100,noise=0.3)
plot(data$x1, data$x2, col=data$target)

#######  Mock Kernel matrix ###### 
rbf.kernel <- function(x, y, sigma=10) {
  return(exp(-sum((x-y)^2)^2/(2*sigma^2)))
}

kernelFun <- "id"
kernel <- switch(kernelFun,
                 rbf=rbfdot(sigma = 1),
                 id=polydot(degree = 1, scale = 1, offset = 0))

K <- kernelMatrix(rbf, data[,1:2])

###### Option A: alternative func rbf & kernel matrix ###### 

## Euclidean Distance  
#1.|| a || = sqrt(aDOTa), 
#2. d(x,y) = || x - y || = sqrt((x-y)DOT(x-y))
#3. aDOTb = sum(a*b)
d<-function(x,y){
  aux=x-y
  dis=sqrt(sum(aux*aux))
  return(dis)
}

##Radial Basis Function Kernel
# 1.K(x,x')=exp(-q||x-x'||^2) where ||x-x'|| is could be defined as the
# euclidian distance and 'q' it's the gamma parameter
rbf<-function(x,y,q=0.2){
  aux<-d(x,y)
  rbfd<-exp(-q*(aux)^2)
  return(rbfd)
}


###### yet another alternative func rbf & kernel matrix ###### 
rbfkernelMx <-function(X, sigma = 1, Y = NULL){
    # test if X is a matrix
    if(!is.matrix(X))
    {
      print("X must be a matrix containing samples in its rows")
      return()
    }
    # test if sigma is a number and > 0
    if(length(sigma) != 1 || sigma <= 0)
    {
      print("sigma must be a number > 0 specifying the rbf-kernel width")
      return()
    }
    if(!is.null(Y))
    {
      # test if Y is a matrix
      if(!is.matrix(Y))
      {
        print("Y must be a matrix containing samples in its rows or NULL if it should not be used")
        return()
      }
      # test if vectors in X and Y have same dimension
      if(ncol(X) != ncol(Y))
      {
        print("The samples in the rows of X and Y must be of same dimension")
        return()
      }
    }
    
    n <- nrow(X) # number of samples in X
    
    if(is.null(Y))
    {
      # calculate distance matrix
      XtX <- tcrossprod(X)
      XX <- matrix(1, n) %*% diag(XtX)
      D <- XX - 2*XtX + t(XX) # distance matrix
    }
    else
    {
      m <- nrow(Y) # number of samples in Y
      # calculate distance matrix (between vectors of X and Y)
      XX <- matrix(apply(X ^ 2, 1, sum), n, m)
      YY <- matrix(apply(Y ^ 2, 1, sum), n, m, byrow = TRUE)
      XY <- tcrossprod(X, Y)
      D <- XX - 2*XY + YY
    }
    
    # calculate rbf-kernel matrix
    K <- exp(-D/(2*sigma))
    
    return(K)
  }


######################################################
##############  KERNEL K-MEANS FUNCTION ##############
######################################################

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

#######  Test functions ###### 

data(spirals)
data <- spirals
rm(spirals)

##### Option A 
#calculating the kernel matrix
kernelmatrix=matrix(0,nrow(data),nrow(data))
for(i in 1:nrow(data)){
  for(j in 1:nrow(data)){
    kernelmatrix[i,j]=rbf(data[i,1:(ncol(data)-1)],data[j,1:(ncol(data)-1)],q=32)
  }
}
K <- kernelmatrix
rm(kernelmatrix)

##### Option B
K <- rbfkernelMx ( X= data, sigma = 0.08)

###### agains normal kmeans with kernel=Identity
clusters <- 2
# Mock assignments
assignments <- matrix(0, nrow = nrow(K), ncol = clusters)
part <- nrow(K)/clusters +1
rassign <- sapply(1:nrow(K)/part, function(x) floor(x)+1)
for(i in 1:nrow(K)) {
  assignments[i,rassign[i]] <- 1
}
# Own Kernel Kmeans
Z <- kernel.kmeans(K, clusters = clusters, iter.max = 10,
                   assignments = assignments)
target <- apply(Z, 1, which.max)
# Regular Kmeans
centers <- (t(assignments) %*% as.matrix(data[,1:2])) / apply(assignments, 2, sum)
kmns <- kmeans(data[,1:2], iter.max = 10,centers = centers)
# See cluster distribution on mock data to see their are equal
#plot(data$x1, data$x2, col=target)
#plot(data$x1, data$x2, col=kmns$cluster)
# plot spiral data 
plot(data[,1], data[,2], col=target, main="Custom Kernel K-means sigma 0.08")
plot(data[,1], data[,2], col=kmns$cluster, main="K-means")

