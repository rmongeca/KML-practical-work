library(caret)
library(mcclust)
library(combinat)

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
    # Reference classes
    ref <- apply(R, 1, which.max) %>% factor(labels=labels)
    # Permutation of columns
    perm <- permn(1:ncol(Z))
    # Get error rate of each permutation and take
    # best permutation to compute performance
    best.error.rate <- Inf
    for(p in perm) {
        Zp <- Z[,p]
        pred <- apply(Zp, 1, which.max) %>% factor(labels=labels)
        # Generate confusion matrix
        cm <- confusionMatrix(pred, ref)
        # Get error rate
        error.rate <- mean(pred != ref)
        if(best.error.rate > error.rate) {
            best.error.rate <- error.rate
            best.cm <- cm
            best.p <- p
            best.pred <- pred
        }
    }
    names(best.error.rate) <- "Error rate"
    # Get variation of information measure
    vi <- vi.dist(best.pred, ref)
    names(vi) <- "Variation of Information"
    return(list(
        confusion.matrix=best.cm$table,
        error.rate=best.error.rate,
        accuracy=best.cm$overall[1],
        kappa=best.cm$overall[2],
        sensitivity=best.cm$byClass[1],
        precision=best.cm$byClass[5],
        recall=best.cm$byClass[6],
        vi=vi,
        permutation=best.p
    ))
}
