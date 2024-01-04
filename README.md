# MSGMM

R software scripts for fitting Gaussian mixture models to multiple data samples simultaneously, using the EM algorithm.
It makes use of C++ code via the R package **Rcpp**.
MSGMM also takes advantage of the open source C++ linear algebra library 'Armadillo' via the R package **RcppArmadillo**.
MSGMM also uses the **data.table** package.
These packages need to be installed first in order to use MSGMM.

The method implemented in MSGMM is detailed in: "Computationally efficient multi-sample flow cytometry data analysis using Gaussian mixture models", *To be submitted.*

Philip Rutten, Tim Robert Mocking, Jacqueline Cloos, Wessel van Wieringen, Costa Bachas

### Running the code

MSGMM is not an R package.
Store the msgmm.R file in the same directory as your R script (or specify the full path name).
The msgmm.R file in turn includes C++ code from external files using the ```sourceCpp("...")``` function.
Make sure that the filename inside sourceCpp() has the right path, and the C++ files are stored in the corresponding location.
Source the R script:

\> source("msgmm.R")

The main funtion ```fitMSGMM()``` fits a GMM to multiple samples, using the EM algorithm, and returns estimated model parameters. 
fitMSGMM iterates data files multiple times and at each iteration opens and closes them sequentially. 
Files should be preprocessed and in CSV format. 
See the reference manual for detailed usage information.

