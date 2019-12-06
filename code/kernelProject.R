#------------------------------- KML-Project -------------------------------#
#PROCEDURE....= 1. Tune Kernel for Kernel-Kmeans and Spectral Clustering    #
#               2. Test performance comparing Variation of Inforation with  #
#                    the true partition of the data                         #
#INPUT DATA...= Pre-processed Reutres data                                  #
#R VERSION....= R version 3.5.1 (2018-07-02) Feather Spray                  # 
#AUTHOR.......= Maria Gkotsopoulou & Ricard Monge Calvo                     #
#CREATED......= December 2019                                               #
#---------------------------------------------------------------------------#

######################################################
#################     LIBRARIES     ##################
######################################################

library(dplyr)
library(kernlab)
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

# We keep the topics acq, trade and ship 
topics <- c("acq", "trade", "ship")
# Own colnames
cols <- c("Topic", "Document")
# Change cols and filter topics
colnames(reutrain) <- cols
colnames(reutest) <- cols
reutrain <- reutrain %>% filter(Topic %in% topics)
reutest <- reutest %>% filter(Topic %in% topics)

# Get topic counts to select topic
topic.counts <- reutrain %>% group_by(Topic) %>% summarise(Count=n()) %>% arrange(desc(Count))
# Get number of instances and topics
N <- dim(reutrain)[1]
k <- length(unique(reutrain$Topic))
# Shuffle data
set.seed(42)
reutrain <- reutrain[sample(1:N, N),]

######################################################
#################   MODEL TUNING    ##################
######################################################
# As string kernel we use Kernlab's version of Normalized Boudrange kernel
# which computes the appearances of subsequences of length less or
# equal than a hyperparameter *length* and then computes the scalar
# product between those appearances vectors.
# To tune the hyperparameter *length* we use 10-fold CV with a range spaces
# going from 2-21 length.
cv <- 10
# Get gridspace to use in parameter tuning
grid.space <- data.frame(length=2:10)
# Get control function to use in parameter tuning by choosing measure in 
# assigment.perfomance function
metric <- "vi"
control <- function(Z, R) {
   perf <- assignment.performance(Z, R)
   mes <- perf[metric]
   return(mes)
}
kernel.tune <- function(
  data, # Dataset to use in the tuning procedure
  cv=10, # Folds to use in K-fold CV
  grid.space, # Dataframe with parameters, where each row is turned into a list
              # of parameters for method function.
  method, # Function to compute Assignment object to compare to Reference.
  reference, # Reference object to use with control function to compare to model
  control, # Control function to compute performance, which recieves the.
          # output of method and the reference object and returns a performance
          # measure that is a number.
  clusters=k # Number of cluster for method
) {
  out.cv <- grid.space %>% mutate(CV=0)
  # List parameters from grid.space
  parameters <- split(grid.space, seq(nrow(grid.space)))
  for(i in 1:nrow(grid.space)) {
    # Get parameters grom grid
    param <- parameters[[i]]
    # Additional parameters for kernel definition
    param$type <- "spectrum"
    param$normalized <- TRUE
    # Get kernel function with parameters
    skernel <- do.call(stringdot, param)
    # Init df to keep cv performances
    aux.cv <- rep(0,cv)
    # Split data into folds
    folds <- createFolds(data,cv)
    for(j in 1:cv) {
      # For each fold create KernelMatrix
      fold.ind <- folds[[j]]
      train <- data[-fold.ind]
      K <- kernelMatrix(skernel, train)
      Z <- method(K, clusters=clusters)
      perf <- control(Z, reference)
      aux.cv[j] <- perf
    }
    out.cv$CV[i] <- mean(aux.cv)
  }
  return(out.cv)
}

