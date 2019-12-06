library(caret)
library(mcclust)

## Assignment vector to matrix
# Utility function to get an assignment matrix from a vector of
# classes for every individual.
assignment.matrix <- function(classes) {
    Z <- matrix(0, nrow = length(classes), ncol = length(unique(classes)))
    for(i in 1:length(classes)) {
        Z[i,classes[i]] <- 1
    }
    return(Z)
}

## Assignment matrix to vector
# Utility function to get an assignment vector from an indicator
# matrix of class assignments.
# Optionally pass arguments for factor function.
assignment.vector <- function(Z, ...) {
    classes <- apply(Z, 1, which.max)
    classes <- factor(classes, ...)
    return(classes)
}

## Assignment performance measure
# Takes an input and reference assignment and computes performance
# measures, and a confusion matrix.
assignment.performance <- function(
    Z, # Assignment to measure as indicator matrix
    R, # Reference assignment as indicator matrix
    labels=1:ncol(R)
) {
    # Number of clusters
    clusters <- ncol(Z)
    pred <- apply(Z, 1, which.max) %>% factor(labels=labels)
    ref <- apply(R, 1, which.max) %>% factor(labels=labels)
    # Generate confusion matrix
    cm <- confusionMatrix(pred, ref)
    # Get error rate
    error.rate <- (cm$table[1,2]+cm$table[2,1])/nrow(R)
    names(error.rate) <- "Error rate"
    # Get variation of information measure
    vi <- vi.dist(pred, ref)
    names(vi) <- "Variation of Information"
    return(list(
        confusion.matrix=cm$table,
        error.rate=error.rate,
        accuracy=cm$overall[1],
        kappa=cm$overall[2],
        sensitivity=cm$byClass[1],
        precision=cm$byClass[5],
        recall=cm$byClass[6],
        vi=vi
    ))
}
