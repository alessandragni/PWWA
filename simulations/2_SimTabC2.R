##### Simulation Study in Web Appendix C (mentioned in Section 5) #####
# Code used for producing Web Table C.2

library(mets)
library(Rcpp)
library(timereg)
library(doMC)

setwd("~/PWWA")

source("simulations/utilsTabC2.R")

RcppArmadillo::armadillo_set_number_of_omp_threads(1)
cc=detectCores()
registerDoMC(cc)

###

set.seed(1)
true10 <- onerun(1, 400000, cens=NULL, time=10, k = 3)
true10


n <- 1000
nsim <- 1000
k = 3


for(cens in c(0.01, 0.03, 0.05)){
  res <- foreach (i=0:nsim) %dopar% { 
    result <- NULL
    while (is.null(result)) {
      try({
        result <- onerun(i, n, time = 10, k = k, cm = ~strata(A, gL), 
                         cens = cens)
      }, silent = TRUE)
    }
    result
  }

  
  write.table(ana(res, true = true10$coef), 
              file=paste("simulations/intermediate_results/TabC1TabC2/C2results_1_",
                         cens,
                         "_insidek",
                         k,
                         ".txt", sep = ""), 
              row.names=c("A=1", "A=0"))
  
}




