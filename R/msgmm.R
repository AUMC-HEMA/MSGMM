library(Rcpp)
library(RcppArmadillo)
library(data.table)

sourceCpp("logstepper_logL_v2.cpp")
sourceCpp("clusterer_v2.cpp")
sourceCpp("loglikelihood.cpp")
sourceCpp("loglikelihoodValues.cpp")

fitMSGMM <- function(samplefiles, 
                     K, 
                     usecols=NULL,
                     init_means=NULL,
                     subsamples=50,
                     init_size=1e4,
                     convergence_threshold=1e-3,
                     max_iter=50,
                     gamma=1, 
                     lambda=0.01, 
                     pooled=FALSE){ 
  
  M <- length(samplefiles)
  
  if (is.null(usecols)){
    p <- ncol(fread(file=samplefiles[1]))
    usecols<-1:p
  } else {
    p <- length(usecols)
  }

  ##############################################################################
  # Initialize parameters

  if (is.null(init_means)){
    
    # Create subsample
    
    Msub <- subsamples
    
    if (length(samplefiles) < Msub){
      Msub <- length(samplefiles)
    }
    
    subsample <- matrix(NA, Msub * init_size, p)
    
    cat("Model",K,"start building subsample from",Msub,"out of",M,"samples.","\n")
    
    for (s in 1:Msub) {
      fcsFile <- samplefiles[s]
      Y1 <- fread(file=fcsFile)
      Y1 <- as.matrix(Y1)[,usecols]
      
      a <- (s - 1) * init_size + 1
      b <- a + init_size - 1
      
      if (nrow(Y1) > init_size){
        subsample[a:b,] <- Y1[sample(nrow(Y1), init_size, replace=FALSE),]
      } else {
        subsample[a:b,] <- Y1[sample(nrow(Y1), init_size, replace=TRUE),]
      }
    }
    
    cat("Subsample size:",Msub * init_size,'\n')

    means <- kmeans(subsample, K, iter.max=1000, algorithm="MacQueen")$centers
    
  } else {
    means <- init_means
  }
  
  ##############################################################################
  
  Sigmas <- t(matrix(rep(diag(p)*100, K), p, K*p))*gamma**2
  
  convergence <- 1
  iter <- 1
  
  ##############################################################################
  # EM loop
  
  if (pooled==TRUE){
    
    weights <- rep(1/K, length(samplefiles), K)
    
    while (convergence > convergence_threshold) {
      cat("Model",K,"iteration",iter,format(Sys.time(),usetz=TRUE),"\n")
      
      # Accumulators
      A0 <- rep(0,K) # matrix(0, 1, K) 
      A1 <- matrix(0, K, p)
      A2 <- matrix(0, K*p, p)
      
      logL <- 0
      
      N <- 0
      
      for (s in 1:M) {
        fcsFile <- samplefiles[s]
        Y1 <- fread(file=fcsFile)
        Y1 <- as.matrix(Y1)[,usecols]
        
        logL <- logL + logstep(Y1, weights, means, Sigmas, A0, A1, A2)
        
        N <- N + nrow(Y1)
      }
      
      if (iter > 1) {
        if (is.na(logL)){
          cat("Error: loglikelihood is",logL,'\n')
          convergence <- 0
        } else {
          convergence <- abs((logL - logL_old) / logL_old)
        }
      }
      
      logL_old <- logL
      iter <- iter + 1
      if (iter > max_iter){
        convergence <- 0
        if (is.na(a)){
          convergence <- 0
        }
      }
      
      # M-step
      weights <- A0 / N
      for (k in 1:K) {
        means[k,] <- A1[k,] / A0[k]
        Sigmas[(k-1)*p + c(1:p),] <- A2[(k-1)*p + c(1:p),] / A0[k] - t(tcrossprod(means[k,], means[k,]))
      }
      
      # Prevent non-invertible matrices
      for (k in 1:K) {
        Sigmas[(k-1)*p + c(1:p),] <-
          (1 - lambda)*Sigmas[(k-1)*p + c(1:p),] + lambda*sum(diag(Sigmas[(k-1)*p + c(1:p),]))/p*diag(p)
      }
    }
    
    d <- list("weights" = weights, "means" = means, "covariances" = Sigmas)
    
    return(d)
    
  } else {
    
    weights <- matrix(1/K, length(samplefiles), K)
    
    while (convergence > convergence_threshold) {
      
      cat("Model",K,"iteration",iter,format(Sys.time(),usetz=TRUE),"\n")
      
      # Accumulators
      A0 <- matrix(0, length(samplefiles), K) # rep(0,K)
      A1 <- matrix(0, K, p)
      A2 <- matrix(0, K*p, p) 
      
      logL <- 0
      
      for (s in 1:M) {
        fcsFile <- samplefiles[s]
        Y1 <- fread(file=fcsFile)
        Y1 <- as.matrix(Y1)[,usecols] 
        
        A0_ <- A0[s,]
        
        logL <- logL + logstep(Y1, weights[s,], means, Sigmas, A0_, A1, A2)
        
        A0[s,] <- A0_
        
        weights[s,] <- A0[s,] / nrow(Y1) 
      }
      
      if (iter > 1) {
        if (is.na(logL)){
          cat("Error: loglikelihood is",logL,'\n')
          convergence <- 0
        } else {
          convergence <- abs((logL - logL_old) / logL_old)
        }
      }
      
      logL_old <- logL
      iter <- iter + 1
      if (iter > max_iter){
        convergence <- 0
      }
      
      # M-step
      for (k in 1:K) {
        means[k,] <- A1[k,] / colSums(A0)[k]
        Sigmas[(k-1)*p + c(1:p),] <- A2[(k-1)*p + c(1:p),] / colSums(A0)[k] - t(tcrossprod(means[k,], means[k,]))
      }
      
      # Prevent non-invertible matrices
      for (k in 1:K) {
        Sigmas[(k-1)*p + c(1:p),] <-
          (1 - lambda)*Sigmas[(k-1)*p + c(1:p),] + lambda*sum(diag(Sigmas[(k-1)*p + c(1:p),]))/p*diag(p)
      }
    }
    
    d <- list("weights" = weights, "means" = means, "covariances" = Sigmas)
    
    return(d)
  }
}

################################################################################

predictLabels <- function(X, weights, means, covariances){
  
  if (is.null(dim(weights))==FALSE){
    if (dim(weights)[1] > 1){
      weights <- colSums(weights)/dim(weights)[1]
    }
  }
  
  dataClusters <- getClusters(X, weights, means, covariances)
  
  return(dataClusters)
}

################################################################################

getLoglikelihood <- function(samplefiles, usecols=NULL, weights, means, covariances){
  
  M <- length(samplefiles)
  
  if (is.null(usecols)){
    p <- ncol(fread(file=samplefiles[1]))
    usecols<-1:p
  } 
  
  logL <- 0
  
  for (s in 1:M) {
    fcsFile <- samplefiles[s]
    cat("Analysing",fcsFile,"\n")
    Y1 <- fread(file=fcsFile)
    Y1 <- as.matrix(Y1)[,usecols]
    
    logL <- logL + getLoglike(Y1, weights, means, covariances)
  }
  
  return(logL)
}

################################################################################

getLoglikelihoodValues <- function(samplefiles, usecols=NULL, weights, means, covariances){
  
  M <- length(samplefiles)
  
  if (is.null(usecols)){
    p <- ncol(fread(file=samplefiles[1]))
    usecols<-1:p
  } 
  
  logLvalues <- c()
  
  for (s in 1:M) {
    fcsFile <- samplefiles[s]
    cat("Analysing",fcsFile,"\n")
    Y1 <- fread(file=fcsFile)
    Y1 <- as.matrix(Y1)[,usecols]
    
    logLvalues <- c(logLvalues,getLoglikeVals(Y1, weights, means, covariances))
  }
  
  return(logLvalues)
}

