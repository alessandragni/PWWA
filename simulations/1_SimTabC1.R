##### Simulation Study in Web Appendix C (mentioned in Section 5) #####
# Code used for producing Web Table C.1


library(timereg)
library(Rcpp)
library(mets)
library(doMC)

setwd("~/PWWA")

# Import functions
Rcpp::sourceCpp("simulations/functionsTabC1.cpp")
source("simulations/utilsTabC1.R")



cc = detectCores()
registerDoMC(cc)


nsim <- 1000 
n = 1000 

t = 10
k = 3


for (S in 1:6){
  for (CENS in c(0.01, 0.03, 0.05)){ # 

    print(S)
    print(CENS)
    
    res10 <- foreach (i = 1:nsim, .errorhandling = "pass") %dopar% onerun(i, n = n, t = t, k = k, S = S, CENS = CENS)
    
    coef <- do.call("rbind", lapply(res10, function(x) x$coef) )
    se.coef <- do.call("rbind", lapply(res10, function(x) x$se.coef) )
    
    write.table(cbind(coef[,1:2], se.coef, coef[,3:6]), 
                file=paste("simulations/intermediate_results/TabC1TabC2/C1results_",
                           S,
                           "_",
                           CENS,
                           "_insidek",
                           k,
                           ".txt", sep = ""), 
                row.names=FALSE, col.names=c('OS10_1', 'OS10_0', 'se10_1', 'se10_0',
                                             'PlugIn1', 'PlugIn0', 'DeBias1', 'DeBias0'))
  }
}

