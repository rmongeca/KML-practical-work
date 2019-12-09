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
library(tidyr)
library(ggplot2)
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
reutrain <- read.delim("data/r52-train-stemmed.txt", header = FALSE,
                       stringsAsFactors = FALSE)
reutest <- read.delim("data/r52-test-stemmed.txt", header = FALSE,
                      stringsAsFactors = FALSE)
reuters <- rbind(reutrain ,reutest) 
rm(reutrain, reutest)
# Own colnames
cols <- c("Topic", "Document")
# Change cols
colnames(reuters) <- cols
# Get topic counts to select topic
topic.counts <- reuters %>% group_by(Topic) %>% summarise(Count=n()) %>% arrange(desc(Count))
# We keep the topics crude, trade and ship 
topics <- c("interest", "crude", "coffee")
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
    Z <- kernel.kmeans(K, clusters=nclust, nrandomSet = nCV)
  }
  if(method=="spectralClustering"){
    Z <- spectralClustering(S=K, clusters=nclust, nrandomSet= nCV )
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
              mutate(type="boundrange",
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
clustMethod <- c("spectralClustering","kernel.kmeans")
gridSpace <- data.frame(length=c(2,5,7,10,12,15,20,25))
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
corpus <- VCorpus(VectorSource(reuters$Document))
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
rm(corpus, dt.mat, dist.matrix)

######################################################
############ Compare &  Best Tune    #################
######################################################

dfplot <- cbind(kernelTuneCV, km=perfKM$vi) %>%
            dplyr::select(-c(type,normalized ))%>%
            tidyr::gather("method","vi" ,CVkk:km)%>%
            mutate(method = ifelse(method=="CVkk", "kernelKmeans", 
                                   ifelse(method=="CVspecc", "spectralClustering",
                                          "kmeans")))

ggplot(dfplot, aes(length,  vi, color = method))+
geom_line() +
geom_point()+
labs(x="string length", y="variation of information") + 
theme_bw() +
theme( panel.background = element_rect(colour = "black", size=0.5),
       panel.grid.minor = element_blank(),
       panel.grid.major = element_blank(),
      axis.text.y = element_text(size=8),
      axis.text.x = element_text(size=8),
      axis.title.y=element_text(size=9),
      axis.title.x=element_text(size=9),
      plot.title = element_text(size=10),
      legend.position="bottom",
      legend.title = element_blank(),
      legend.text = element_text( size=7),
      legend.key.height=unit(0.5,"line"),
      legend.key.width=unit(0.5,"line"))

######################################################
############ Run Methods   #################
######################################################

nrandomS <- 10
## best lambda
lambda_kk <- dfplot %>% filter(method == "kernelKmeans") %>% 
                slice(which.min(vi)) %>% pull(length)
## compute Kernel Matrix with best lambda
sk_kk <- stringdot(length= lambda_kk,type = "spectrum", normalized = TRUE)
K <- kernelMatrix(sk_kk, tuneData)
## do clustering method using customs functions
kernkk <- kernel.kmeans(K, clusters=k, nrandomSet = nrandomS)

## best lambda
lambda_spec <- dfplot %>% filter(method == "spectralClustering") %>% 
                slice(which.min(vi)) %>% pull(length)
## compute Kernel Matrix with best lambda
sk_spec <- stringdot(length= lambda_spec,type = "spectrum", normalized = TRUE)
K <- kernelMatrix(sk_spec, tuneData)
## do clustering method using customs functions
spec <- spectralClustering(S=K, clusters=k, nrandomSet= nrandomS )

######################################################
############ Performance &  Compare  #################
######################################################

perfKernkk <- assignment.performance(kernkk,ref)

perfSpec <- assignment.performance(spec,ref)

perfDF <- tibble(kernelkmeans = perfKernkk$vi,specc = perfSpec$vi,kmeans = perfKM$vi)%>%
            add_row(kernelkmeans=perfKernkk$error.rate,specc=perfSpec$error.rate, kmeans=perfKM$error.rate) %>%
            add_row(kernelkmeans=perfKernkk$accuracy,specc=perfSpec$accuracy, kmeans=perfKM$accuracy) %>%
            add_row(kernelkmeans=perfKernkk$sensitivity,specc=perfSpec$sensitivity, kmeans=perfKM$sensitivity) %>%
            add_row(kernelkmeans=perfKernkk$precision,specc=perfSpec$precision, kmeans=perfKM$precision) %>%
            add_row(kernelkmeans=perfKernkk$recall,specc=perfSpec$recall, kmeans=perfKM$recall) 


rownames(perfDF) <- c("vi", "error.rate", "accuracy", "sensitivity","precision", "recall")


saveRDS(perfDF, file = "performanceMeasures")



