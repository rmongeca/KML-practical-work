# KML-practical-work
Practical work for the subject Kernel-based Machine Learning of the UPC MIRI master course.

## Library dependencies
Install the following packages:
* dplyr, kernlab, tm, proxy, tidyr, ggplot2, caret, mcclust, combinat

## Run instructions
In an R Studio session open the Rproject KML-practical-work. The main code is code/kernelProject, in which 
the other 3 codes are sourced. Specifically, the kernelKmeans includes our implementation of the *kernel k-means* method while, the SpectralClustering includes our implementation of the *spectral clustering* method. Finally, AssignmentPerformance includes utility function to compute performance measures used to perform comparison across methods and with respect to the true partition. 

Furthermore, for convenience purposes we attach the dataset used to run all the experiments that can be found in the data folder. The author created a train and test partition that we source separately but later merge. The code kernelProject includes the following sections:

* data load
* model tuning Functions
* model tuning run
* k-means process including preprocessing (dtm)
* tuning comparison
* model run using the best hyperparameter
* performance comparison