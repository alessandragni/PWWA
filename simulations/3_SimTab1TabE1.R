##### Simulation Study in Section 5 #####
# Code used for producing Table 1 and Web Table E.1

# Import libraries
library(mets)
library(doMC)
library(RcppArmadillo)

setwd("~/PWWA")

# Import functions
source("simulations/utils.R")

# Load data and set up simulation parameters
load("data/hhosp.rda")

dfactor(hhosp) <- treat~trt_ab

cox1 <- phreg(Event(start,stop,status==1)~factor(treat)+cluster(id),data=hhosp)
coxd <- phreg(Event(start,stop,status==3)~factor(treat)+cluster(id),data=hhosp)

Lam12 <- cox1$cumhaz
LamD <- coxd$cumhaz

Lam11 <- cbind(c(0,1,4.31),c(0,2.5,3.45))
Lam13 <- cbind(c(0,1,4.31),c(0,0.5,3.45))
Lams <- list(Lam11,Lam12,Lam13)

set.seed(1)
seeds <- sample(1:10^6, size = 5001, replace = FALSE)

##### Run some simulations #####

RcppArmadillo::armadillo_set_number_of_omp_threads(1)
cc=detectCores()
registerDoMC(cc)

######  Tables 1 and E.1 ######
beta11 <- -0.3
betad1 <- -0.3
time = 3

# Parse input parameters
args <- commandArgs(trailingOnly = TRUE)
varz <- as.numeric(args[1])
cens <- as.numeric(args[2])

type = 2
dep = 1
scale1 = 1
scaled = 1

n <- 1000
nsim <- 5000 

ntrue <- 10000
nsimtrue <- 500


# set values
# depvals <- c(1) # dependence: 1 = shared frailty (v=1), 4 = frailty only recurrent events (v=0)
# typevals <- c(2)
# cens_vals <- c(1/4,2/4)
# theta_vals <- c(0.5,1,2)
# scale1_vals <- c(1)
# scaled_vals <- c(1)

# for (type in typevals) 
#   for (dep in depvals) 
#     for (varz in theta_vals) 
#       for (cens in cens_vals) 
#         for (scale1 in scale1_vals) 
#           for (scaled in scaled_vals) {
            
            print(c(varz,cens))
            
            # computation of the "mean","sd","mean-se",("power") 
            reslZ <- resl <- list()
            outtotZ <- outtot <- c()
            
            # a model that included $\widetilde{N}(r-)$, $L$, 
            # and $Z$ for censoring augmentation and both $L$ and $Z$ for the mean ratio model.
            # for Table 1
            
            reslZ <- foreach (i=1:nsim) %dopar% onerunN(i,n,var.z=varz,time=time,dep=dep,
                                                        scale1=scale1,scaled=scaled,
                                                        cens=cens,trans=0.33,beta11=beta11,betad1=betad1,
                                                        Lam1=Lams[[type]],LamD=LamD,seeds=seeds)
            
            # models that included only $\widetilde{N}(r-)$ and 
            # $L$ for censoring augmentation and only $L$ for the mean ratio model
            # for Table E.1
            
            resl <- foreach (i=1:nsim) %dopar% onerunNS(i,n,var.z=varz,time=time,dep=dep,
                                                        scale1=scale1,scaled=scaled,
                                                        cens=cens,trans=0.33,beta11=beta11,betad1=betad1,
                                                        Lam1=Lams[[type]],LamD=LamD,seeds=seeds)

            mmZ <- ana(reslZ)
            outtotZ <- cbind(n,cens,dep,varz,scale1,mmZ$summary)
            
            mm <- ana(resl)
            outtot <- cbind(n,cens,dep,varz,scale1,mm$summary)
            
            
            # computation of the true values 
            # as the mean of simulations to obtain also the "bias", "mse" and the "coverage"
            
            restrueZ <- restrue <- list()
            
              
            # for Table 1
            restrueZ <- foreach (i=1:nsimtrue) %dopar% onerunN(i,ntrue,var.z=varz,time=time,dep=dep,
                                                               scale1=scale1,scaled=scaled,
                                                               cens=NULL,trans=0.33,beta11=beta11,betad1=betad1,
                                                               Lam1=Lams[[type]],LamD=LamD,seeds=seeds)
            # for Table E.1
            restrue <- foreach (i=1:nsimtrue) %dopar% onerunNS(i,ntrue,var.z=varz,time=time,dep=dep,
                                                               scale1=scale1,scaled=scaled,
                                                               cens=NULL,trans=0.33,beta11=beta11,betad1=betad1,
                                                               Lam1=Lams[[type]],LamD=LamD,seeds=seeds)
            
            
            # save results
            outtab1 <- list(resl=reslZ,outtot=outtotZ,true=restrueZ)
            outtabE1 <- list(resl=resl,outtot=outtot,true=restrue)
            
            save(outtab1, file=paste0("simulations/intermediate_results/Tab1TabE1/Tab1_varz",varz,"_cens",cens,".rda")) 
            save(outtabE1, file=paste0("simulations/intermediate_results/Tab1TabE1/TabE1_varz",varz,"_cens",cens,".rda"))
            
#          }
