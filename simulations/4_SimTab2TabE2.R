##### Simulation Study in Section 5 #####
# Code used for producing Table 2 and Web Table E.2

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

######  Table 2 and Web Table E.2 ######

time = 3

# Parse input parameters
args <- commandArgs(trailingOnly = TRUE)
type <- as.numeric(args[1])
dep <- as.numeric(args[2])
varz <- as.numeric(args[3])
scaled <- as.numeric(args[4])

cens = 1/4
scale1 = 1

n <- 1000
nsim <- 5000

ntrue <- 10000
nsimtrue <- 100


# set values
# depvals <- c(1,4) # dependence: 1 = shared frailty (v=1), 4 = frailty only recurrent events (v=0)
# typevals <- c(2,3,1)
# theta_vals <- c(1,2)
# scaled_vals <- c(1,4)
# cens_vals <- c(1/4)
# scale1_vals <- c(1)

# for (type in typevals) 
#   for (dep in depvals) 
#     for (varz in theta_vals) 
#       for (cens in cens_vals) 
#         for (scale1 in scale1_vals) 
#           for (scaled in scaled_vals) {
            
            print(c(type,dep,varz,scaled))
            
            # computation of the "mean","sd","mean-se",("power") 
            resl <- reslnull <- list()
            outtot <- outtotnull <- c()
              
            # for Table 2
            resl <- foreach (i=1:nsim) %dopar% onerunNR(i,n,var.z=varz,time=time,dep=dep,
                                                        scale1=scale1,scaled=scaled,
                                                        cens=cens,trans=0.33,beta11=-0.3,betad1=-0.3,
                                                        Lam1=Lams[[type]],LamD=LamD,seeds=seeds)
            
            # for Table E.2
            reslnull <- foreach (i=1:nsim) %dopar% onerunNR(i,n,var.z=varz,time=time,dep=dep,
                                                            scale1=scale1,scaled=scaled,
                                                            cens=cens,trans=0.33,beta11=0,betad1=0,
                                                            Lam1=Lams[[type]],LamD=LamD,seeds=seeds)

            mm <- anaR(resl)
            outtot <- cbind(n,type,dep,varz,scaled,mm$summary)
            
            mmnull <- anaR(reslnull)
            outtotnull <- cbind(n,type,dep,varz,scaled,mmnull$summary)
            
            
            # computation of the true values 
            # as the mean of simulations to obtain also the "bias", "mse" and the "coverage"
            restrue <- reslnulltrue <- list()
              
            # for Table 2
            restrue <- foreach (i=1:nsimtrue) %dopar% onerunNR(i,ntrue,var.z=varz,time=time,dep=dep,
                                                               scale1=scale1,scaled=scaled,
                                                               cens=NULL,trans=0.33,beta11=-0.3,betad1=-0.3,
                                                               Lam1=Lams[[type]],LamD=LamD,seeds=seeds)
            
            # for Table E.2
            reslnulltrue <- foreach (i=1:nsimtrue) %dopar% onerunNR(i,ntrue,var.z=varz,time=time,dep=dep,
                                                                    scale1=scale1,scaled=scaled,
                                                                    cens=NULL,trans=0.33,beta11=0,betad1=0,
                                                                    Lam1=Lams[[type]],LamD=LamD,seeds=seeds)
            
            # save results
            
            outtab2 <- list(resl=resl,outtot=outtot,true=restrue)
            outtabE2 <- list(resl=reslnull,outtot=outtotnull,true=reslnulltrue)
            
            save(outtab2, file=paste0("simulations/intermediate_results/Tab2TabE2/Tab2_type",type,"_dep",
                                      dep,"_varz",varz,"_scaled",scaled,".rda")) 
            save(outtabE2, file=paste0("simulations/intermediate_results/Tab2TabE2/TabE2_type",type,"_dep",
                                       dep,"_varz",varz,"_scaled",scaled,".rda"))
            
#          }
