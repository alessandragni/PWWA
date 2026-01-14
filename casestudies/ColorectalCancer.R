#####  Case Study based on Colorectal cancer data - Web Appendix D #####
# Code used for producing Figures D.1 and D.2 in Web Appendix D

###### Import package ######
library(mets) 

###### Load data ######
load("data/colorectal.rda")

###### Preprocessing ######
rr <- colorectal
rrid <- countID(rr)
rr <- cbind(rr,rrid)
rr$treatment <- factor(rr$treatment)
rr <- count.history(rr,status="new.lesions",types=1)
rr$treattime <- (rr$lbnr__id==1)*1
###
rr$status <- rr$new.lesions
rr <- dtransform(rr,status=2,state==1)
dtable(rr,~status)
rr$rid <- rr$reverseCountid
rr$trt <- rr$treatment


####### Web Appendix D - Figure D.1 (i) #######
pdf("casestudies/results/FigD1i-MM_treatment.pdf", width=7, height=5)
xr <- phreg(Surv(time0,time1,new.lesions)~strata(treatment)+survival::cluster(id),data=rr)
dr <- phreg(Surv(time0,time1,state)~strata(treatment)+survival::cluster(id),data=rr)
out <- recurrentMarginal(xr,dr)
plot(out,se=TRUE,ylab="Marginal Mean",xlab='Time (years)',col=1:2, lwd=3)
stratnames <- paste('treatment', unique(rr$treatment), sep = ":")
legend("topleft", legend=stratnames, col=1:2, lwd=3, lty=1:2)
dev.off()

###### Web Appendix D - Figure D.1 (ii) ######
pdf("casestudies/results/FigD1ii-MM_prevresection.pdf", width=7, height=5)
xr <- phreg(Surv(time0,time1,new.lesions)~strata(prev.resection)+survival::cluster(id),data=rr)
dr <- phreg(Surv(time0,time1,state)~strata(prev.resection)+survival::cluster(id),data=rr)
out <- recurrentMarginal(xr,dr)
plot(out,se=TRUE,ylab="Marginal Mean",xlab='Time (years)',col=3:4, lwd=3)
stratnames <- paste('prev.resection', unique(rr$prev.resection), sep = ":")
legend("topleft", legend=stratnames, col=3:4, lwd=3, lty=1:2)
dev.off()

###### Web Appendix D - Figure D.1 (iii)  ######
pdf("casestudies/results/FigD1iii-MM_Age.pdf", width=7, height=5)
xr <- phreg(Surv(time0,time1,new.lesions)~strata(age)+survival::cluster(id),data=rr)
dr <- phreg(Surv(time0,time1,state)~strata(age)+survival::cluster(id),data=rr)
out <- recurrentMarginal(xr,dr)
plot(out,se=TRUE,ylab="Marginal Mean",xlab='Time (years)', col=1:3, lwd=3, lty=c(1,2,6))
stratnames <- paste('age', c('<60 years', '60-69 years', '>69 years'), sep = ":")
legend("topleft", legend=stratnames, col=1:3, lwd=3, lty=c(1,2,6))
dev.off()

###### Web Appendix D - Figure D.1 (iv) ######
pdf("casestudies/results/FigD1iv-MM_whoPS.pdf", width=7, height=5)
xr <- phreg(Surv(time0,time1,new.lesions)~strata(who.PS)+survival::cluster(id),data=rr)
dr <- phreg(Surv(time0,time1,state)~strata(who.PS)+survival::cluster(id),data=rr)
out <- recurrentMarginal(xr,dr)
plot(out,se=TRUE,ylab="Marginal Mean",xlab='Time (years)',col=6:8, lwd=3, lty=c(1,2,6))
stratnames <- paste('who.PS', unique(rr$who.PS), sep = ":")
legend("topleft", legend=stratnames, col=6:8, lwd=3, lty=c(1,2,6))
dev.off()




###### Computation of the estimators of interest ######
res <- resR <- c()
for (tt in seq(0.5,2.9,by=0.1)) {
  dd0 <- WA_recurrent(Event(time0,time1,status)~treatment+cluster(id),
                      rr,time=tt,death.code=2,
                      augmentR=~age+who.PS+prev.resection,
                      augmentC=~age+who.PS+prev.resection) #, trans = 1/3)
  ee <- estimate(coef=dd0$ET$riskDR$riskDR,vcov=dd0$ET$riskDR$var.riskDR)$coefmat
  eeR <- dd0$RAW$ratio.means$coefmat
  tdif <- abs(diff(ee[,1])/sum(ee[,2]^2)^.5)
  pval <- 2*(1-pnorm(tdif))
  tdifR <- abs(diff(eeR[,1])/sum(eeR[,2]^2)^.5)
  pvalR <- 2*(1-pnorm(tdifR))
  res <- rbind(res,c(tt,ee[,1],ee[1,3:4],ee[2,3:4],pval))
  resR <- rbind(resR,c(tt,eeR[,1],eeR[1,3:4],eeR[2,3:4],pvalR))
}


###### Plots ######
plotres <- function(res, ylab) {
  matplot(res[,1],res[,2:3],type="l",lwd=3,ylim=c(0.2,1.2),xlab="Time (years)",ylab=ylab)
  plotConfRegion(res[,1],res[,4:5],col=1)
  plotConfRegion(res[,1],res[,6:7],col=2)
  if (ncol(res)==8) {
    sigp <- (res[,8]<0.05)
    points(res[sigp,1],res[sigp,2],pch="*",cex = 2.5, col=1)
    points(res[sigp,1],res[sigp,3],pch="*",cex = 2.5, col=2)
  }
  legend("bottomright",c("treatment:S","treatment:C"),lty=1:2,col=1:2,lwd=2.5)
}

####### Web Appendix D - Figure D.2 (i) #######
pdf("casestudies/results/FigD2i-colorectal-wwa-PWWA.pdf", width=5.5, height=5)
plotres(res, ylab = "PWWA estimand")
dev.off()

####### Web Appendix D - Figure D.2 (ii) #######
pdf("casestudies/results/FigD2ii-colorectal-wwa-EWWA.pdf", width=5.5, height=5)
plotres(resR, ylab = "EWWA estimand")
dev.off()

