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
#################     DATA LOAD     ##################
######################################################
## We are going to use a slightly-processed version of the famous Reuters news
## articles dataset. All articles with no Topic annotations are dropped. The
## text of each article is converted to lowercase, whitespace is normalized to
## single-spaces. Only the first term from the Topic annotation list is retained
## (some articles have several topics assigned).  

## The resulting dataset is a list of pairs (Topic, News).
## The news text is the column "Content" and its category is the column "Topic".
## The goal is test clustering methods and evaluate the variation of information
## of their partitions (assignments) with the true partition of the data.

## Note that we can directly read the compressed version (reuters.txt.gz). 
# There is no need to unpack the gz file; for local files R handles unpacking.
reuters <- read.table("data/reuters.txt.gz", header=T)

# R originally loads this as factor, so needs fixing
reuters$Content <- as.character(reuters$Content)

length(reuters$Topic)

table(reuters$Topic)

N <- dim(reuters)[1]

set.seed(42)
reuters <- reuters[sample(1:N, N),]