library(MASS)
library(ellipse)

################################################################################
# Demonstration of MSGMM analysis on simulated sample files

setwd("/Users/philip/Documents/postdoc/vumc/code/Rcpp/test")
source("msgmm.R")

################################################################################
# First, generate four bivariate sample files with four clusters.
# One sample is much stronger expressed in one of the clusters. 
# The means and covariances are fixed, except the mixing proportions.
# We store the samples in separate CSV files.

setwd("/Users/philip/Documents/postdoc/vumc/code/Rcpp/demo-analysis")

K <- 4 # number of clusters
p <- 2 # number of features
S <- 4 # number of samples

pimatrix <- matrix(c(0.3, 0.3, 0.3, 0.2,
                     0.3, 0.3, 0.3, 0.2,
                     0.3, 0.3, 0.3, 0.2,
                     0.1, 0.1, 0.1, 0.4),4,4)

means <- matrix(c(1,1,-1,-1,1,-1,-1,1),4,2)

covariances <- matrix(c(1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1)*0.1,8,2)

for (s in 1:S) {
  clusterlabels <- sample(K, 1000, replace = TRUE, pimatrix[s,])
  cells <- matrix(NA, 1000, 2)
  for (k in unique(clusterlabels)) {
    nk <- which(clusterlabels==k)
    cells[nk,] <- mvrnorm(length(nk),
                          means[k,],
                          covariances[(k-1)*p + c(1:p),])
  }
  
  write.csv(cells,paste("sim-",toString(s),".csv",sep = ''),row.names = FALSE)
}

################################################################################
# Now, analyse the stored sample files.
patients <- c("sim-1.csv","sim-2.csv","sim-3.csv","sim-4.csv")

# Fit the MSGMM to the samples
# Note that setting gamma << 1 is important here!
model_params <- fitMSGMM(patients, 
                         K = 4, 
                         init_size = 1e3,
                         convergence_threshold=1e-4, 
                         gamma = 0.1)

# Create a scatter plot with fitted model parameters.
Y <- read.csv(patients[1])

plot(Y[,1], Y[,2], pch=19, cex=1, col="grey")

for (k in 1:K){
  points(model_params$means[k,1],
         model_params$means[k,2], pch=19, cex=1, col="black")
  text(model_params$means[k,1],model_params$means[k,2], 
       toString(k), cex=2, 0, col="black")
  lines(ellipse(model_params$covariances[(k-1)*p + 1,2],
                centre=model_params$means[k,],
                level=0.5), lwd=2, lty=1, col="black")
}

# Create a clustered heatmap of the pi-matrix
# You see the one sample with aberrant cell expressions 
heatmap(model_params$weights,labRow = patients)
