getSubsample <- function(files, K, usecols, init.files, init.size, seed){

  if(!is.null(seed)){
    set.seed(seed)
  }

  if (length(files) < init.files){
    cat("init.files is smaller than input files, setting equal to input files \n")
    init.files <- length(files)
  }
  subsample <- matrix(NA, init.files * init.size, length(usecols))
  cat("Model",K,"start building subsample from",init.files,"out of",length(files),"samples.","\n")
  for (s in 1:init.files) {
    Y1 <- data.table::fread(file=files[s])
    Y1 <- as.matrix(Y1)[,usecols]
    a <- (s - 1) * init.size + 1
    b <- a + init.size - 1
    if (nrow(Y1) > init.size){
      subsample[a:b,] <- Y1[sample(nrow(Y1), init.size, replace=FALSE),]
    } else {
      subsample[a:b,] <- Y1[sample(nrow(Y1), init.size, replace=TRUE),]
    }
  }
  cat("Subsample size:",init.files * init.size,'\n')
  return(subsample)
}


#' Fit multi-sample Gaussian Mixture Model (MSGMM)
#' 
#' @param files Character vector with CSV filenames
#' @param K Integer specifying number of Gaussian components
#' @param usecols Vector specifying the variables (columns) to use
#' @param init.means Initial values of means. (K-means used if not specified)
#' @param init.files Number of files to use for means initialization
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
MSGMM <- function(files, K, usecols = NULL, init.means = NULL, init.files = 50,
                  init.size = 1e4, seed = NULL, tol = 1e-3, max.iter = 50,
                  gamma = 1, lambda = 0.01, pooled = FALSE){ 

  if (is.null(usecols)){
    p <- ncol(data.table::fread(file=files[1]))
    usecols<-1:p
  } else {
    p <- length(usecols)
  }

  # Initialize parameters
  if (is.null(init.means)){
    subsample <- getSubsample(files, K, usecols, init.files, init.size, seed)
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
  if (pooled==TRUE){
    while (convergence > tol && iter <= max.iter) {
      cat("Model",K,"iteration",iter,format(Sys.time(),usetz=TRUE),"\n")
      # Initialize accumulators
      A0 <- rep(0,K) # matrix(0, 1, K) 
      A1 <- matrix(0, K, p)
      A2 <- matrix(0, K*p, p)
      logL <- 0
      N <- 0
      
      # E-step
      for (s in 1:length(files)) {
        Y1 <- data.table::fread(file=files[s])
        Y1 <- as.matrix(Y1)[,usecols]
        logL <- logL + logstep(Y1, weights, means, covariances, A0, A1, A2)
        N <- N + nrow(Y1)
      }
      # Assess convergence
      if (iter > 1) {
        if (is.na(logL)){
          stop("Log-likelihood is ", logL)
        } else {
          convergence <- abs((logL - logL_old) / logL_old)
        }
      }
      logL_old <- logL
      iter <- iter + 1
      # M-step
      weights <- A0 / N
      for (k in 1:K) {
        means[k,] <- A1[k,] / A0[k]
        covariances[(k-1)*p + c(1:p),] <- A2[(k-1)*p + c(1:p),] / A0[k] - t(tcrossprod(means[k,], means[k,]))
      }
      # Prevent non-invertible matrices
      for (k in 1:K) {
        covariances[(k-1)*p + c(1:p),] <-
          (1 - lambda)*covariances[(k-1)*p + c(1:p),] + lambda*sum(diag(covariances[(k-1)*p + c(1:p),]))/p*diag(p)
      }
    }
  } else {
    
    while (convergence > tol && iter <= max.iter) {
      cat("Model",K,"iteration",iter,format(Sys.time(),usetz=TRUE),"\n")
      # Initialize accumulators
      A0 <- matrix(0, length(files), K) # rep(0,K)
      A1 <- matrix(0, K, p)
      A2 <- matrix(0, K*p, p) 
      logL <- 0
      
      # E-step
      for (s in 1:length(files)) {
        Y1 <- data.table::fread(file=files[s])
        Y1 <- as.matrix(Y1)[,usecols] 
        A0_ <- A0[s,]
        logL <- logL + logstep(Y1, weights[s,], means, covariances, A0_, A1, A2)
        A0[s,] <- A0_
        weights[s,] <- A0[s,] / nrow(Y1) 
      }
      # Assess convergence
      if (iter > 1) {
        if (is.na(logL)){
          stop("Log-likelihood is ", logL)
        } else {
          convergence <- abs((logL - logL_old) / logL_old)
        }
      }
      logL_old <- logL
      iter <- iter + 1
      # M-step
      for (k in 1:K) {
        means[k,] <- A1[k,] / colSums(A0)[k]
        covariances[(k-1)*p + c(1:p),] <- A2[(k-1)*p + c(1:p),] / colSums(A0)[k] - t(tcrossprod(means[k,], means[k,]))
      }
      # Prevent non-invertible matrices
      for (k in 1:K) {
        covariances[(k-1)*p + c(1:p),] <-
          (1 - lambda)*covariances[(k-1)*p + c(1:p),] + lambda*sum(diag(covariances[(k-1)*p + c(1:p),]))/p*diag(p)
      }
    }
  }
  output <- list("weights" = weights, "means" = means, "covariances" = covariances)
  class(output) <- "MSGMM"
  return(output)
}


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


getLoglikelihood <- function(files, usecols = NULL, params){
  if (is.null(usecols)){
    p <- ncol(data.table::fread(file=files[1]))
    usecols<-1:p
  } 
  logL <- 0
  for (s in 1:length(files)) {
    cat("Analysing",files[s],"\n")
    Y1 <- data.table::fread(file=files[s])
    Y1 <- as.matrix(Y1)[,usecols]
    logL <- logL + getLoglike(Y1, params$weights, params$means, params$covariances)
  }
  return(logL)
}


getLoglikelihoodValues <- function(files, usecols = NULL, params){
  if (is.null(usecols)){
    p <- ncol(data.table::fread(file=files[1]))
    usecols<-1:p
  } 
  logLvalues <- c()
  for (s in 1:length(files)) {
    cat("Analysing",files[s],"\n")
    Y1 <- data.table::fread(file=files[s])
    Y1 <- as.matrix(Y1)[,usecols]
    logLvalues <- c(logLvalues,getLoglikeVals(Y1, params$weights, 
                                              params$means, params$covariances))
  }
  return(logLvalues)
}
