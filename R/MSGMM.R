#' Wrapper for reading files from different sources (CSV, FCS)
#' 
#' @param file File path
#' @param usecols Vector specifying the variables (columns) to use
#'
#' @return Numeric matrix
readData <- function(file, usecols){
  if (substring(file, nchar(file) - 3, nchar(file)) == ".csv") {
    mat <- as.matrix(data.table::fread(file = file, check.names = FALSE,
                                       select = usecols))
  }
  if (substring(file, nchar(file) - 3, nchar(file)) == ".fcs") {
    # Read only relevant columns (if character)
    if (is.character(usecols)){
      patt <- paste(usecols, collapse = "|")
      patt <- paste0("^(", patt, ")$")
      mat <- flowCore::read.FCS(file, column.pattern = patt,
                                truncate_max_range = FALSE)@exprs
    } else {
      mat <- flowCore::read.FCS(file, truncate_max_range = FALSE)@exprs[, usecols]
    }
  }
  return(mat)
}


#' Subsample data from a list of CSV or FCS files
#' 
#' @param init.files Character vector CSV or FCS filenames to use for K-means
#' @param usecols Vector specifying the variables (columns) to use
#' @param init.size Number of data points to sample from each file for means initialization 
#' @param seed Random seed for subsampling during means initialization 
#'
#' @return Matrix with subsampled data
subsampleFiles <- function(init.files, usecols, init.size, seed){
  # Set seed for random sampling
  if(!is.null(seed)){
    set.seed(seed)
  }
  matList <- list()
  for (s in seq_along(init.files)) {
    X <- readData(init.files[s], usecols)
    rows <- X[sample(nrow(X), min(nrow(X), init.size), replace = FALSE), ]
    matList[[s]] <- rows
  }
  subsample <- do.call(rbind, matList)
  colnames(subsample) <- NULL
  return(subsample)
}


#' Fit multi-sample Gaussian Mixture Model (MSGMM)
#' 
#' @param files Character vector with CSV or FCS filenames
#' @param K Integer specifying number of Gaussian components
#' @param usecols Numeric or character vector specifying the columns to use
#' @param init.means Initial values of means. (K-means used if not specified)
#' @param init.files Character vector with CSV or FCS filenames to use for K-means
#' @param init.size Number of data points to sample from each file for means initialization 
#' @param seed Random seed for subsampling during means initialization 
#' @param tol Convergence threshold for stopping criterion 
#' @param max.iter Maximum number of EM iterations
#' @param gamma Starting value of component variance matrices 
#' @param lambda Regularization parameter for component covariance estimation 
#' @param pooled Logical flag to enable/disable multi-sample EM algorithm
#'
#' @return List of model parameters ("weights", "means", and "covariances")
#'
#' @export
MSGMM <- function(files, K, usecols, init.means = NULL, init.files = NULL,
                  init.size = 1e4, seed = NULL, tol = 1e-3, max.iter = 50,
                  gamma = 1, lambda = 0.01, pooled = FALSE){
  
  # Initialize parameters
  p <- length(usecols)
  if (is.null(init.means)){
    if (is.null(init.files)){
      init.files <- files
    }
    subsample <- subsampleFiles(init.files, usecols, init.size, seed)
    means <- stats::kmeans(subsample, K, iter.max=1000, algorithm="MacQueen")$centers
  } else {
    means <- init.means
  }
  if (pooled == TRUE){
    weights <- rep(1/K, length(files), K)
  } else {
    weights <- matrix(1/K, length(files), K)
  }
  covariances <- t(matrix(rep(diag(p)*100, K), p, K*p))*gamma**2
  
  # EM loop
  convergence <- 1
  iter <- 1
  while (convergence > tol && iter <= max.iter) {
    cat("Model",K,"iteration",iter,format(Sys.time(),usetz=TRUE),"\n")
    # Initialize accumulators
    logL <- 0
    N <- 0
    if (pooled){
      A0 <- rep(0,K) # matrix(0, 1, K) 
      A1 <- matrix(0, K, p)
      A2 <- matrix(0, K*p, p)
    } else {
      A0 <- matrix(0, length(files), K) # rep(0,K)
      A1 <- matrix(0, K, p)
      A2 <- matrix(0, K*p, p) 
    }
    # E-step
    for (s in 1:length(files)) {
      X <- readData(files[s], usecols)
      N <- N + nrow(X)
      if (pooled){
        logL <- logL + logstep(X, weights, means, covariances, A0, A1, A2)
      } else {
        A0_ <- A0[s,]
        logL <- logL + logstep(X, weights[s,], means, covariances, A0_, A1, A2)
        A0[s,] <- A0_
        weights[s,] <- A0[s,] / nrow(X) 
      }
    }
    # Assess convergence
    if (iter > 1) {
      if (is.na(logL)){
        stop("Log-likelihood is ", logL)
      } else {
        convergence <- abs((logL - logL_old) / logL_old)
      }
    }
    # Update parameters (M-step)
    logL_old <- logL
    iter <- iter + 1
    if (pooled){
      weights <- A0 / N
      for (k in 1:K) {
        means[k,] <- A1[k,] / A0[k]
        covariances[(k-1)*p + c(1:p),] <- A2[(k-1)*p + c(1:p),] / A0[k] - t(tcrossprod(means[k,], means[k,]))
      }
    } else {
      for (k in 1:K) {
        means[k,] <- A1[k,] / colSums(A0)[k]
        covariances[(k-1)*p + c(1:p),] <- A2[(k-1)*p + c(1:p),] / colSums(A0)[k] - t(tcrossprod(means[k,], means[k,]))
      }
    }
    # Prevent non-invertible matrices
    for (k in 1:K) {
      covariances[(k-1)*p + c(1:p),] <-
        (1 - lambda)*covariances[(k-1)*p + c(1:p),] + lambda*sum(diag(covariances[(k-1)*p + c(1:p),]))/p*diag(p)
    }
  }
  metadata <- list("files" = files, "K" = K, "usecols" = usecols, 
                   "init.means" = init.means, "init.files" = init.files,
                   "init.size" = init.size, "seed" = seed, "tol" = tol, 
                   "max.iter" = max.iter, "gamma" = gamma, "lambda" = lambda, 
                   "pooled" = pooled, "n" = N, "p" = p)
  output <- list("weights" = weights, "means" = means, "covariances" = covariances,
                 "metadata" = metadata)
  class(output) <- "MSGMM"
  return(output)
}


#' Assign cluster labels to new data
#' 
#' @param X Numeric matrix
#' @param params MSGMM model parameter output
#'
#' @return Cluster labels
#' 
#' @export
predictLabels <- function(X, params){
  weights <- params$weights
  if (is.null(dim(weights))==FALSE){
    if (dim(weights)[1] > 1){
      weights <- colSums(weights)/dim(weights)[1]
    }
  }
  dataClusters <- getClusters(X, weights, params$means, params$covariances)
  return(dataClusters)
}


#' Calculate log-likelihood
#' 
#' @param files Character vector with CSV or FCS filenames
#' @param params MSGMM model parameter output
#'
#' @return Log-likelihood values
#' 
#' @seealso \code{\link{getLoglikelihoodValues}}
#'
#' @export
getLoglikelihood <- function(files, params){
  usecols <- params$metadata$usecols
  logL <- 0
  for (s in 1:length(files)) {
    cat("Analysing",files[s],"\n")
    X <- readData(files[s], usecols)
    if (params$metadata$pooled){
      logL <- logL + getLoglike(X, params$weights, params$means, params$covariances)
    } else {
      if (sum(basename(model_params$metadata$files) == basename(files)) != length(model_params$metadata$files)){
        message("Warning: input files different from files used in fit!\nSample-specific weights are invalid, defaulting to GMM...")
        # Average the model weights
        weights <- colMeans(model_params$weights)
        logL <- logL + getLoglike(X, weights, params$means, params$covariances)
      } else {
        logL <- logL + getLoglike(X, params$weights[s,], params$means, params$covariances)
      }
    }
  }
  return(logL)
}


#' Calculate log-likelihood values
#' 
#' @param files Character vector with CSV or FCS filenames
#' @param params MSGMM model parameter output
#'
#' @return Vector with log-likelihood values
#' 
#' @seealso \code{\link{getLoglikelihood}}
#'
#' @export
getLoglikelihoodValues <- function(files, params){
  usecols <- params$metadata$usecols
  logLvalues <- c()
  for (s in 1:length(files)) {
    cat("Analysing",files[s],"\n")
    X <- readData(files[s], usecols)
    logLvalues <- c(logLvalues,getLoglikeVals(X, params$weights, 
                                              params$means, params$covariances))
  }
  return(logLvalues)
}


#' Calculate the Bayesian Information Criterion (BIC)
#' 
#' @param logL Log-likelihood value of fitted model
#' @param params MSGMM model parameter output
#'
#' @return BIC score
#' 
#' @seealso \code{\link{getLoglikelihood}}
#'
#' @export
calculateBIC <- function(logL, params){
  n <- params$metadata$n
  K <- params$metadata$K
  p <- params$metadata$p
  s <- length(params$metadata$files)
  pooled <- params$metadata$pooled
  if (pooled){
    nParams <- (K * p) + K * (p * (p + 1) / 2) + (K - 1)
  } else {
    nParams <- (K * p) + K * (p * (p + 1) / 2) + (K - 1) * s
  }
  BIC <- -2 * logL + nParams * log(n)
  return(BIC)
}


#' Get features from each file based on model parameters
#' 
#' @param files Character vector with CSV or FCS filenames
#' @param params MSGMM model parameter output
#' @param type Type of features to extract ("counts" or "percentages")
#'
#' @return Matrix with features per file
#' 
#' @export
getFeatures <- function(files, params, type = "percentages"){
  usecols <- params$metadata$usecols
  output <- matrix(nrow = length(files), ncol = nrow(params$means))
  for (s in seq(1, length(files))){
    cat("Analysing",files[s],"\n")
    X <- readData(files[s], usecols)
    labels <- predictLabels(X, params)
    features <- table(factor(labels, levels = 1:nrow(params$means)))
    if (type == "percentages"){
      features <- features / sum(features)
    }
    output[s,] <- features
  }
  rownames(output) <- files
  return(output)
}
