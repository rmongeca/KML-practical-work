#------------------------------- KML-Project -------------------------------#
#PROCEDURE....= 1. Tune Kernel for Kernel-Kmeans and Spectral Clustering    #
#               2. Test performance comparing Variation of Inforation with  #
#                    the true partition of the data                         #
#INPUT DATA...= Pre-processed Reuters data                                  #
#R VERSION....= R version 3.5.1 (2018-07-02) Feather Spray                  # 
#AUTHOR.......= Maria Gkotsopoulou & Ricard Monge Calvo                     #
#CREATED......= December 2019                                               #
#---------------------------------------------------------------------------#

######################################################
#################     LIBRARIES     ##################
######################################################

library(dplyr)
library(kernlab)
library(tm)
library(proxy)
set.seed(42)

######################################################
#################   LOAD FUNCTIONS  ##################
######################################################

source("code/AssignmentPerformance.R")
source("code/kernelKmeans.R")
source("code/SpectralClustering.R")


######################################################
#################     DATA LOAD     ##################
######################################################
# Text files containing one document per line.
# Each document is composed by its class and its terms.
# Each document is represented by a "word" representing the document's class,
# a TAB character and then a sequence of "words" delimited by spaces,
# representing the terms contained in the document.
# Train/test set
reutrain <- read.delim("data/r52-train-no-stop.txt", header = FALSE,
                       stringsAsFactors = FALSE)
reutest <- read.delim("data/r52-test-no-stop.txt", header = FALSE,
                      stringsAsFactors = FALSE)
reuters <- rbind(reutrain ,reutest) 
rm(reutrain, reutest)
# Get topic counts to select topic
topic.counts <- reuters %>% group_by(Topic) %>% summarise(Count=n()) %>% arrange(desc(Count))
# We keep the topics crude, trade and ship 
topics <- c("crude", "trade", "ship")
# Own colnames
cols <- c("Topic", "Document")
# Change cols and filter topics
colnames(reuters) <- cols
reuters <- reuters %>% filter(Topic %in% topics) %>% 
              mutate(Topic=factor(Topic))
# Get number of instances and topics
N <- dim(reuters)[1]
k <- length(unique(reuters$Topic))

######################################################
#########   MODEL TUNING Functions   #################
######################################################
# As string kernel we use Kernlab's version of Normalized Boudrange kernel
# which computes the appearances of subsequences of length less or
# equal than a hyperparameter *length* and then computes the scalar
# product between those appearances vectors.
# To tune the hyperparameter *length* we use 10-fold CV with a range spaces
# going from 2-21 length.

# Get control function to use in parameter tuning by choosing measure in 
# assigment.perfomance function
metric <- "vi"
# Control function to compute performance, which recieves the.
# output of method and the reference object and returns a performance
# measure that is a number.
control <- function(Z, R) {
  perf <- assignment.performance(Z, R)
  mes <- unlist(perf[metric])
  return(mes)
}

clusteringMethod <- function(
  method,
  K, # kernelMatrix
  nclust, # Number of cluster for method
  nCV){
  if(method=="kernel.kmeans"){
    Z <- kernel.kmeans(K, clusters=nclust,nrandomSet = nCV)
  }
  if(method=="spectralClustering"){
    Z <- spectralClustering(S=K, clusters=nclust,nrandomSet= nCV )
  }
  return(Z)
}

kernel.tune <- function(
  data, # Dataset to use in the tuning procedure
  clusters ,
  cv, # Folds to use in K-fold CV
  grid.space, # Dataframe with parameters, where each row is turned into a list
              # of parameters for method function.
  # Function to compute Assignment object to compare to Reference.
  clMethod= c("kernel.kmeans", "spectralClustering"), 
  reference # Reference object to use with control function to compare to model
) {
  out.cv <- grid.space %>% 
              mutate(type="spectrum",
                     normalized = TRUE,
                     CVkk = 0, 
                     CVspecc = 0)
  # List parameters from grid.space
  parameters <- split(grid.space, seq(nrow(grid.space)))
  for(i in 1:nrow(grid.space)) {
    # Get parameters grom grid
    param <- parameters[[i]]
    # Get kernel function with parameters
    skernel <- do.call(stringdot, param)
    # Init df to keep cv performances
    aux.cv <- rep(0,cv)
    K <- kernelMatrix(skernel, data)
    for(x in 1:length(clMethod)){
      Z <- clusteringMethod(method=clMethod[x], K , nclust =clusters,nCV=cv )
      perf <- control(Z, reference)
      if(clMethod[x]=="spectralClustering"){
        out.cv$CVspecc[i] <- perf
      }
      if(clMethod[x]=="kernel.kmeans"){
        out.cv$CVkk[i] <- perf 
      }
    }
  }
  return(out.cv)
}


######################################################
##################   MODEL TUNING    #################
######################################################
tuneData <- reuters$Document
clustMethod <- c("spectralClustering")
gridSpace <- data.frame(length=c(2,5,10,15,20))
ref <- assignment.matrix(reuters$Topic)

kernelTuneCV <- kernel.tune(data = tuneData ,
                            clusters = k,
                            cv=10,
                            grid.space =  gridSpace,
                            clMethod= clustMethod,
                            reference =ref )

######################################################
##################   K Means    #################
######################################################

# Transform data to DocumentTermMatrix for regular k-means
# use the tf–idf(term frequency–inverse document frequency) 
# instead of the frequencies of the term as entries, 
# tf-idf measures the relative importance of a word to a document.
# and reduce dimensions allowing only up to 0.95 matrix sparsity
corpus <- VCorpus(VectorSource(reuters$Document))%>%
            tm::tm_map(stemDocument)
corpus <- tm_map(corpus, PlainTextDocument)
dt.mat <- DocumentTermMatrix(corpus ,
                             control = list(weighting = weightTfIdf)) %>% 
            removeSparseTerms(sparse = 0.95) 
dist.matrix = proxy::dist(as.matrix(dt.mat), method = "cosine")
km <- kmeans(dist.matrix,centers=k,nstart = 10, iter.max =100)

classes <- km$cluster
assignmentsKM <- matrix(0, nrow = length(classes), ncol = k)
for(i in 1:length(classes)) {
  assignmentsKM[i,classes[i]] <- 1
}
perfKM <- assignment.performance(assignmentsKM, ref)


## plot 

## best lambdas

## do clustering method using customs functions

## with Z do assignment.performance
                            
  

