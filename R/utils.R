#' Splits the concatenated covariance matrix into a list ordered by component
splitCovariances <- function(params){
  K <- nrow(params$means)
  nRows <- nrow(params$covariances)
  rowsPerK <- nRows / K
  covList <- lapply(1:K, function(i) {
    start <- (i - 1) * rowsPerK + 1
    end <- i * rowsPerK
    params$covariances[start:end, ]
  })
  return(covList)
}


#' Calculate Kullback–Leibler (KL) divergence between two components
#' 
#' @param mu1 Vector of means per dimension for first component
#' @param covMat1 Covariance matrix of first component
#' @param mu2 Vector of means per dimension for second component
#' @param covMat2 Covariance matrix of second component
#'
#' @return Vector with log-likelihood values
calculateKL <- function(mu1, covMat1, mu2, covMat2) {
  # Ensure covMat1 and covMat2 are symmetric
  covMat1 <- (covMat1 + t(covMat1)) / 2
  covMat2 <- (covMat2 + t(covMat2)) / 2
  # Use true inverse if covMat1 is invertible, else generalized inverse
  if (det(covMat1) > 1e-10) {
    covMat1inv <- solve(covMat1)
  } else {
    covMat1inv <- MASS::ginv(covMat1)
  }
  traceTerm <- sum(diag(covMat1inv %*% covMat2))
  meanTerm <- t(mu2 - mu1) %*% covMat1inv %*% (mu2 - mu1)
  K <- length(mu1)
  logTerm <- log(det(covMat2) / det(covMat1))
  KL <- 0.5 * (traceTerm + meanTerm - K + logTerm)
  return(KL)
}


#' Calculate Jeffreys divergence between two components
calculateJeffrey <- function(mu1, covMat1, mu2, covMat2) {
  kl1 <- calculateKL(mu1, covMat1, mu2, covMat2)
  kl2 <- calculateKL(mu2, covMat2, mu1, covMat1)
  return(kl1 + kl2)
}


#' Calculate distance matrix between components based on Jeffreys divergence
#' 
#' @param params MSGMM model parameter output
#' 
#' @return Distance matrix
#'
#' @export
compDist <- function(params){
  K <- nrow(params$means)
  # Split the covariance matrix by component
  covList <- splitCovariances(params)
  distMat <- matrix(0, nrow = K, ncol = K)
  for (k1 in seq(1, K)){
    for (k2 in seq(1, K)){
      mu1 <- params$means[k1,]
      mu2 <- params$means[k2,]
      covMat1 <- covList[[k1]]
      covMat2 <- covList[[k2]]
      # Calculate symmetric KL divergence
      distMat[k1, k2] <- calculateJeffrey(mu1, covMat1, mu2, covMat2)
    }
  }
  return(distMat)
}


#' Get component ellipsoids based on model parameters
#' 
#' @param params MSGMM model parameter output
#' @param dim1 First dimension
#' @param dim2 Second dimension
#' @param conf Confidence level
#' 
#' @return Ellipsoid dimensions (x, y)
#'
#' @export
getEllipsoids <- function(params, dim1, dim2, conf = 0.95){
  covList <- splitCovariances(params)
  K <- nrow(params$means)
  ellipsoids <- list()
  ellipsoids <- vector("list", K)  # Initialize list with predefined length
  for (i in seq_len(K)) {  
    cov2D <- covList[[i]][c(dim1, dim2), c(dim1, dim2)]
    mu2D <- params$means[i, c(dim1, dim2)]
    el <- ellipse::ellipse(cov2D, centre = mu2D, level = conf)
    ellipsoids[[i]] <- el
  }
  return(ellipsoids)
}
