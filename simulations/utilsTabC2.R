
#### Useful functions ####

expit = function(x){ return( exp(x) / (exp(x) + 1)) }


#### Data Simulation ####

sim.data2 = function(sample.size = 1000, CENS = 0.005){# {{{
  n = sample.size
  
  id = 1:n
  
  L = runif(n, min=0, max=1) # covariate
  
  A = rbinom(n, 1, expit(-0.5 + L)) # treatment c(rep(1, n/2), rep(0, n/2))
  
  beta1 = log(2)
  gam1.A = 1
  
  lambda1 = 0.04 * exp(beta1 * L + gam1.A * A)
  lambda2 = 0.02 * exp(beta1 * L + gam1.A * A)
  lambda12 = 0.05 * exp(beta1 * L + gam1.A * A)
  
  T1minT2 = rexp(n, (lambda1 + lambda2))
  delta1 = rbinom(n, size=1, prob = lambda1/(lambda1 + lambda2))
  T2 = numeric(n)
  T2[delta1==0] = T1minT2[delta1==0]
  
  U = rexp(n, lambda12)
  T2[delta1==1] = T1minT2[delta1==1] + U[delta1==1]
  
  theta = 1
  
  if (!is.null(CENS)) {
    lambdaC = CENS * exp(theta * A + 1 * I(L>0.5))
    C = rexp(n, lambdaC)
  } else C <- max(T2)+1
  
  
  Ttilde2 = pmin(T2, C)
  delta2 = (Ttilde2 == T2) * 1
  Ttilde1 = pmin(T1minT2, Ttilde2)
  delta1 = (T1minT2 < Ttilde2) * 1
  
  censbeforeany = ((delta1 == 0 ) & (delta2 == 0))     # 1 whether "censored censored before any event" or not
  censafter1before2 = ((delta1 == 1 ) & (delta2 == 0)) # 1 whether "non-terminal event then censored prior to terminal event" or not
  cens = (censbeforeany | censafter1before2) * 1       # 1 whether censored (at any point) or not
  
  d <- data.frame(
    id = id,
    L = L,
    A = A,
    C = C,
    Ttilde1 = Ttilde1,
    Ttilde2 = Ttilde2,
    delta1 = delta1,
    delta2 = delta2,
    cens = cens
  )
  return(d)
}
# }}}




#### One run ####

dynCensAug <- function(formC,data,augmentC=~+1,response="Yipcw",time=NULL,Z=NULL) { ## {{{
  if (is.null(time)) stop("must give time of response \n")
  data$Y__ <- data[,response]
  varsC <- c("Y__",attr(terms(augmentC), "term.labels"))
  formCC <- update(formC, reformulate(c(".", varsC)))
  cr2 <- phreg(formCC, data = data, no.opt = TRUE, no.var = 1,Z=Z)
  xx <- cr2$cox.prep
  icoxS0 <- rep(0,length(cr2$S0))
  icoxS0[cr2$S0>0] <- 1/cr2$S0[cr2$S0>0]
  S0i <- rep(0, length(xx$strata))
  S0i[xx$jumps + 1] <- icoxS0
  km <- TRUE
  if (!km) {
    cumhazD <- c(cumsumstratasum(S0i, xx$strata, xx$nstrata)$lagsum)
    St <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1 - S0i), xx$strata, xx$nstrata)$lagsum))
  
  nterms <- cr2$p-1
  dhessian <- cr2$hessianttime
  dhessian <-  .Call("XXMatFULL",dhessian,cr2$p,PACKAGE="mets")$XXf
  ###  matrix(apply(dhessian,2,sum),3,3)
  timeb <- which(cr2$cumhaz[,1]<time)
  ### take relevant \sum H_i(s,t) (e_i - \bar e)
  covts <- dhessian[timeb,1+1:nterms,drop=FALSE]
  ### construct relevant \sum (e_i - \bar e)^2
  Pt <- dhessian[timeb,-c((1:(nterms+1)),(1:(nterms))*(nterms+1)+1),drop=FALSE]
  ###  matrix(apply(dhessian[,c(5,6,8,9)],2,sum),2,2)
  gammatt <- .Call("CubeVec",Pt,covts,1,PACKAGE="mets")$XXbeta
  S0 <- cr2$S0[timeb]
  gammatt[is.na(gammatt)] <- 0
  gammatt[gammatt==Inf] <- 0
  Gctb <- St[cr2$cox.prep$jumps+1][timeb]
  augmentt.times <- apply(gammatt*cr2$U[timeb,1+1:nterms,drop=FALSE],1,sum)
  augment.times <- sum(augmentt.times)
  if (!is.null(Z)) {
    Zj <- cr2$cox.prep$Z[cr2$cox.prep$jumps+1,][timeb]
    Xaugment.times <- apply( augmentt.times*Zj,2,sum)
  }
  
  p <- 1
  #### iid magic  for censoring augmentation martingale #
  ### int_0^infty gamma(t) (e_i - ebar(s)) 1/G_c(s) dM_i^c
  xx <- cr2$cox.prep
  nid <- max(xx$id)+1
  jumpsC <- xx$jumps+1
  rr0 <- xx$sign
  S0i <- rep(0,length(xx$strata))
  S0i[jumpsC] <- c(1/(icoxS0*St[jumpsC]))
  S0i[jumpsC] <- icoxS0
  xxtime <- 1*c(xx$time<time)
  
  pXXA <- ncol(cr2$E)-1
  EA <- cr2$E[timeb,-1,drop=FALSE]
  gammasEs <- .Call("CubeMattime",gammatt,EA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
  gammasE <- matrix(0,length(xx$strata),1)
  gammattt  <-    matrix(0,length(xx$strata),pXXA*1)
  jumpsCt <- jumpsC[timeb]
  gammasE[jumpsCt,] <- gammasEs
  gammattt[jumpsCt,] <- gammatt
  gammaEsdLam0 <- apply(gammasE*S0i*xxtime,2,cumsumstrata,xx$strata,xx$nstrata)
  gammadLam0 <-   apply(gammattt*S0i*xxtime,2,cumsumstrata,xx$strata,xx$nstrata)
  XgammadLam0 <- .Call("CubeMattime",gammadLam0,xx$X[,-1,drop=FALSE],pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
  Ut <- Et <- matrix(0,length(xx$strata),1)
  Ut[jumpsCt,] <- augmentt.times
  MGCt <- Ut[,drop=FALSE]-(XgammadLam0-gammaEsdLam0)*c(rr0)
  MGCiid <- apply(MGCt,2,sumstrata,xx$id,nid)
  MGCiid <- MGCiid/nid
  
  nid <- max(cr2$id)+1
  res <- list(MGCiid=MGCiid,gammat=gammatt,augment=augment.times,
              id=cr2$name.id,n=nid)
} ## }}}


onerun <- function(i,n,time=1,cm=~1,k=3,cens=0.005,response="EpT",dcr=1,...) {# {{{
  
  nid <- n
  if (i%%200==0) print(i)
  ###
  dd <- sim.data2(n,CENS=cens)
  dtable(dd,~delta2)
  dd$time <- dd$Ttilde2
  dd$status <- dd$delta2*2
  dcut(dd,breaks=2) <- gL~L
  dd$gL <- (dd$L>0.5)
  ddt <- event.split(dd,cuts="Ttilde1")
  ddt <- dtransform(ddt,status=1,time<Ttilde2)
  ddt <- count.history(ddt,types=1)
  head(ddt)
  ddt$stop <- ddt$time
  rr <- ddt
  rr$statusD <- rr$status
  rr$x <- rr$L
  dfactor(rr) <- treatment~A
  ###
  times <- time
  iddata <- rr
  ## cut after time of interest to compute duration
  iddata <- event.split(iddata,status="statusD",cuts=times+0.1,time="stop",name.start="start")
  rrc <- subset(iddata,stop<=times+0.1)
  rrc$revnrc <- revcumsumstrata(rep(1,nrow(rrc)),rrc$id-1,n)
  ###
  rrc$lasttime <- rrc$stop[rrc$revnrc==1][rrc$id]
  rrc$EpT <- (rrc$Count1/rrc$lasttime)^(1/k) 
  head(rrc)
  ## use this response
  rrcl <- subset(rrc,revnrc==1)
  rrcl$event <- (rrcl$status!=0)*1
  head(rrcl$EpT)
  ###
  ###
  newdata <- rrcl
  
  outae <- binregATE(Event(stop,event)~treatment+x,newdata,cause=1,time=time,treat.model=treatment~x,
                     Ydirect=newdata$EpT,outcome="rmst",model="lin",cens.model=cm) ## ,control=list(stepsize=0.5))
  
  treat.model=treatment~x
  
  treats <- function(treatvar) {# {{{ treatvar <- droplevels(treatvar)
    nlev <- nlevels(treatvar)
    nlevs <- levels(treatvar)
    ntreatvar <- as.numeric(treatvar)
    return(list(nlev = nlev, nlevs = nlevs, ntreatvar = ntreatvar))
  } # }}}
  
  fittreat <- function(treat.model, data, id, ntreatvar, nlev) {# {{{
    if (nlev == 2) {
      treat.model <- drop.specials(treat.model, "cluster")
      treat <- glm(treat.model, data, family = "binomial")
      iidalpha <- lava::iid(treat, id = id)
      lpa <- treat$linear.predictors
      pal <- expit(lpa)
      pal <- cbind(1 - pal, pal)
      ppp <- (pal/pal[, 1])
      spp <- 1/pal[, 1]
    }
    else {
      treat.modelid <- update.formula(treat.model, . ~
                                        . + cluster(id__))
      treat <- mlogit(treat.modelid, data)
      iidalpha <- lava::iid(treat)
      pal <- predictmlogit(treat, data, se = 0, response = FALSE)
      ppp <- (pal/pal[, 1])
      spp <- 1/pal[, 1]
    }
    Xtreat <- model.matrix(treat.model, data)
    tvg2 <- 1 * (ntreatvar >= 2)
    pA <- c(mdi(pal, 1:length(ntreatvar), ntreatvar))
    pppy <- c(mdi(ppp, 1:length(ntreatvar), ntreatvar))
    Dppy <- (spp * tvg2 - pppy)
    Dp <- c()
    for (i in seq(nlev - 1)) Dp <- cbind(Dp, Xtreat * ppp[,
                                                          i + 1] * Dppy/spp^2)
    DPai <- -1 * Dp/pA^2
    out <- list(iidalpha = iidalpha, pA = pA, Dp = Dp, pal = pal,
                ppp = ppp, spp = spp, id = id, DPai = DPai)
    return(out)
  } # }}}
  
  expit <- function(x) 1/(1 + exp(-x))
  idW <- newdata$id
  
  treatsvar <- newdata[,"treatment"]
  treats <- treats(treatsvar)
  
  fitt <- fittreat(treat.model, newdata, idW, treats$ntreatvar, treats$nlev)
  iidalpha0 <- fitt$iidalpha
  wPA <- c(fitt$pA)
  
  newdata$wPA <- wPA[newdata$id]
  
  
  rrc$Yp <- rrc$EpT/wPA[rrc$id]
  rrc$Count1p <- rrc$Count1/wPA[rrc$id]
  rrc$xp <- rrc$x/wPA[rrc$id]
  
  if (dcr==1) {
    dc0 <- dynCensAug(Surv(start,stop,statusD==0)~+cluster(id),subset(rrc,A==0),augmentC=~Count1+x,
                      response="Yp",time=time)
    dc1 <- dynCensAug(Surv(start,stop,statusD==0)~+cluster(id),subset(rrc,A==1),augmentC=~Count1+x,
                      response="Yp",time=time)
  } else {
    dc0 <- dynCensAug(Surv(start,stop,statusD==0)~+cluster(id),subset(rrc,A==0),augmentC=~Count1p+xp,
                      response="Yp",time=time)
    dc1 <- dynCensAug(Surv(start,stop,statusD==0)~+cluster(id),subset(rrc,A==1),augmentC=~Count1p+xp,
                      response="Yp",time=time)
  }
  nn <- table(rrc$A)
  
  
  MGC0 <- sumstrata(dc0$MGCiid,dc0$id-1,nid)*dc0$n
  MGC1 <- sumstrata(dc1$MGCiid,dc1$id-1,nid)*dc1$n
  MGC <- cbind(MGC0,MGC1)
  ccaugment <- apply(MGC,2,sum)
  riskDRC <- outae$riskDR+ccaugment/nid
  
  iidDRC <- outae$riskDR.iid + MGC/nid
  varDRC <- crossprod(iidDRC)
  se.riskDRC <- diag(varDRC)^.5
  
  outaef <- binregATE(Event(stop,event)~treatment*x,newdata,cause=1,time=time,treat.model=treatment~x,
                      Ydirect=newdata$EpT,outcome="rmst",model="lin",cens.model=cm) ## ,control=list(stepsize=0.5))
  summary(outaef)
  
  riskDRCF <- outaef$riskDR+ccaugment/nid
  iidDRCF <- outaef$riskDR.iid+MGC/nid
  varDRCF <- crossprod(iidDRCF)
  se.riskDRCF <- diag(varDRCF)^.5
  
  newdata$wPA <- wPA
  
  out0 <- resmeanIPCW(Event(stop,event)~+1,newdata,cause=1, time=time,type="I",
                      Ydirect=(newdata$A==0)*newdata$EpT/newdata$wPA,model="lin",cens.model=cm) ## ,control=list(stepsize=0.5))
  summary(out0)
  
  out1 <- resmeanIPCW(Event(stop,event)~+1,newdata,cause=1, time=time,type="I",
                      Ydirect=(newdata$A==1)*newdata$EpT/newdata$wPA,model="lin",cens.model=cm) ## ,control=list(stepsize=0.5))
  summary(out1)
  
  Y0 <- (newdata$A==0)*newdata$EpT/newdata$wPA
  Y1 <- (newdata$A==1)*newdata$EpT/newdata$wPA
  mean(Y0)
  mean(Y1)
  out1$cens.weights
  
  
  coef=c(outae$riskDR,outaef$riskDR, riskDRC, riskDRCF,   out0$coef, out1$coef);
  se.coef=c(outae$se.riskDR,outaef$se.riskDR, se.riskDRC,se.riskDRCF, out0$se.coef, out1$se.coef)
  
  names(coef) <- c(rep(c("DR-A+X","DR-A*X","C-DR-A+X","C-DR-A*X"),each=2), rep("Strat-ipcw-A",each=2))
  names(coef) <- c(rep(c("DR-A+X","DR-A*X","C-DR-A+X","C-DR-A*X"),each=2), rep("Strat-ipcw-A",each=2))
  names(se.coef) <- names(coef)
  ###
  res <- list(coef=coef, se.coef=se.coef)
  return(res) 
} # }}}



#### Results Analysis ####

ana <- function(res,true=NULL) {# {{{
  res <- Filter(function(x) !inherits(x, "try-error"), res)
  
  coef <-do.call("rbind",lapply(res,function(x) x$coef) )
  scoef <-do.call("rbind",lapply(res,function(x) x$se.coef) )
  m <- cbind( apply(coef,2,mean), apply(coef,2,mean)-true, apply(coef,2,sd), apply(scoef,2,mean) )
  colnames(m) <- c("mean", "bias", "sd","mean-se") #,"var-ratio")
  
  if (!is.null(true)) {
    covv <- apply(t(t(coef - 1.96*scoef) < true & t(coef + 1.96*scoef) > true),2,mean,na.rm=TRUE)
    m <- cbind(m,covv)
    colnames(m) <- c("mean", "bias", "sd","mean-se","coverage") #,"var-ratio")
  }
  
  # Interleave even rows with odd rows
  odd_indices <- seq(1, nrow(m), by = 2)
  even_indices <- seq(2, nrow(m), by = 2)
  interleaved_indices <- as.vector(rbind(even_indices, odd_indices))
  #interleaved_indices <- interleaved_indices[!is.na(interleaved_indices)] # Remove NA values if any
  
  m <- m[interleaved_indices, ]
  m2 = cbind( m[7:8, 1:5], m[3:4, 1:3])
  colnames(m2) <- c("hat psi", "bias", "sd","mean-se","coverage", "tilde psi", "bias", "sd") 
  
  # apply(round(m2, 3), 1, function(row) {
  #   cat(paste(row, collapse = " & "), "\\\\\n")
  # })
    
  return(m2)
} ## }}}



latex_body <- function(m2, alpha, d = 3) {
  f <- function(x) formatC(x, format = "f", digits = d)
  c(
    sprintf("& \\multirow{2}{*}{$\\alpha = %.2f$} & $A=1$ & %s \\\\",
            alpha, paste(f(unlist(m2[1,])), collapse = " & ")),
    sprintf("& & $A=0$ & %s \\\\",
            paste(f(unlist(m2[2,])), collapse = " & ")),
    "\\cmidrule{2-11}"
  )
}






