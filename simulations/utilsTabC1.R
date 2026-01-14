#### Useful functions ####

expit = function(x){
  return( exp(x) / (exp(x) + 1))
}


g = function(x){
  return(x^(1/k))
}




##### One-step Estimator EFF ##### 


at = function(t, k, df,
              tau12, dLam12,
              tau1, dLam1,
              Lam1.tau1, 
              Lam2.tau1, 
              Lam1.tauc, 
              Lam2.tauc, 
              Lam12.tau1, Lam12.tau12,
              Lam12.tauc, 
              tauc, KC.tauc, dLamC) {
  
  CM1 = CM1_cpp(t, k, df$Ttilde2, df$Ttilde1, df$delta1, 
                tau12, dLam12, 
                tau1, dLam1,
                Lam1.tau1,
                Lam2.tau1,
                Lam1.tauc,
                Lam2.tauc,
                Lam12.tau1, Lam12.tau12,
                tauc, KC.tauc, dLamC)
  
  CM2 = CM2_cpp(t, k, df$Ttilde2, df$Ttilde1, df$delta1, 
                tau12, dLam12, Lam12.tauc, Lam12.tau12,
                tauc, KC.tauc, dLamC)
  
  KCTtilde2 = NULL
  for(i in 1:dim(df)[1]){
    KCTtilde2 = c(KCTtilde2, look_up_f_cpp(df$Ttilde2[i], tauc, KC.tauc[ ,i]))
  }
  
  at = (1 * ((df$Ttilde1 <= t) & (df$delta1 == 1)) * df$delta2 ) / ( g(pmin(df$Ttilde2, t)) * KCTtilde2 ) +  CM1 + CM2
  
  return( at )
}


OneStep = function(t, k, df, pA1,
                   tau12, dLam12.a0, dLam12.a1,
                   tau1, dLam1.a0, dLam1.a1,
                   Lam1.a0.tau1, Lam1.a1.tau1, 
                   Lam2.a0.tau1, Lam2.a1.tau1,
                   Lam1.a0.tauc, Lam1.a1.tauc, 
                   Lam2.a0.tauc, Lam2.a1.tauc,
                   Lam12.a0.tau1, Lam12.a1.tau1,
                   Lam12.a0.tau12, Lam12.a1.tau12, 
                   Lam12.a0.tauc, Lam12.a1.tauc,
                   tauc, KC.a0.tauc, KC.a1.tauc, dLamC.a0, dLamC.a1) {
  
  Ht1 = Heff_cpp(t, k, Lam1.a1.tau1, Lam2.a1.tau1, tau1, dLam1.a1, tau12, dLam12.a1, Lam12.a1.tau1, Lam12.a1.tau12)
  Ht0 = Heff_cpp(t, k, Lam1.a0.tau1, Lam2.a0.tau1, tau1, dLam1.a0, tau12, dLam12.a0, Lam12.a0.tau1, Lam12.a0.tau12)
  
  at1 = at(t, k, df,
           tau12, dLam12.a1,
           tau1, dLam1.a1,
           Lam1.a1.tau1, 
           Lam2.a1.tau1,
           Lam1.a1.tauc, 
           Lam2.a1.tauc,
           Lam12.a1.tau1, Lam12.a1.tau12, Lam12.a1.tauc, 
           tauc, KC.a1.tauc, dLamC.a1)
  
  at0 = at(t, k, df,
           tau12, dLam12.a0, 
           tau1, dLam1.a0, 
           Lam1.a0.tau1, 
           Lam2.a0.tau1, 
           Lam1.a0.tauc, 
           Lam2.a0.tauc, 
           Lam12.a0.tau1, Lam12.a0.tau12, Lam12.a0.tauc, 
           tauc, KC.a0.tauc, dLamC.a0)
  
  debiasing_term1 = (1*(df$A == 1) / pA1) * ( at1 - Ht1 )
  debiasing_term0 = (1*(df$A == 0) / (1-pA1)) * ( at0 - Ht0 ) 
  
  PlugIn1 = mean(Ht1)
  PlugIn0 = mean(Ht0)
  
  OS1 = PlugIn1 + mean(debiasing_term1)
  OS0 = PlugIn0 + mean(debiasing_term0)
  
  se1 = sqrt(mean((Ht1 + debiasing_term1 - OS1)^2)) / sqrt(dim(df)[1])
  se0 = sqrt(mean((Ht0 + debiasing_term0 - OS0)^2)) / sqrt(dim(df)[1])
  
  coef = c(OS1, OS0, PlugIn1, PlugIn0, mean(debiasing_term1), mean(debiasing_term0))
  se.coef = c(se1, se0)
  
  res <- list(coef=coef, se.coef=se.coef)
  return( res )
}




#### Data Simulation ####

sim.data = function(sample.size = 1000, CENS){
  n = sample.size
  
  id = 1:n
  
  L = runif(n, min=0, max=1) # covariate
  #L = 1 * I(L > 0.5)
  
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
  
  lambdaC = CENS * exp(theta * A + 1 * I(L>0.5))
  C = rexp(n, lambdaC)
  
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


#### Data Generation ####

est.f = function(df, t, k, S){
  
  with(df, {
    
    n = dim(df)[1]
    A = df$A
    L = df$L
    tau = df$tau[1]
    
    #### ESTIMATION COX ####
    # Here we will need to fit three different models for dLambda1, dLambda2 and dLambda12
    
    if (S == 3 | S == 4 | S == 6){
      fit.surv1 = cox.aalen(Surv(Ttilde1, delta1) ~ #prop(L) + prop(A) + prop(I(A*L)), 
                              prop(A),
                            max.clust = n, data = df)  
      k1 = dim(fit.surv1$cum)[1]       
      
      g1.a0 = exp( fit.surv1$gamma[1] * rep(0, n))# + fit.surv1$gamma[1] * L + fit.surv1$gamma[3] * I(0*L))
      g1.a1 = exp( fit.surv1$gamma[1] * rep(1, n))# + fit.surv1$gamma[1] * L + fit.surv1$gamma[3] * I(1*L))
    } else {
      fit.surv1 = cox.aalen(Surv(Ttilde1, delta1) ~ prop(L) + prop(A),  
                            max.clust = n, data = df)  
      k1 = dim(fit.surv1$cum)[1]        # extracting the number of time at which event happened!
      
      g1.a0 = exp( fit.surv1$gamma[2] * 0 + fit.surv1$gamma[1] * L )
      g1.a1 = exp( fit.surv1$gamma[2] * 1 + fit.surv1$gamma[1] * L )
    }
    
    
    fit.surv2 = cox.aalen(Surv(Ttilde1, c(1-delta1)*delta2) ~ prop(L) + prop(A), 
                          max.clust = n, data = df) 
    k2 = dim(fit.surv2$cum)[1]       
    
    g2.a0 = exp( fit.surv2$gamma[2] * 0 + fit.surv2$gamma[1] * L)
    g2.a1 = exp( fit.surv2$gamma[2] * 1 + fit.surv2$gamma[1] * L)
    
    
    
    if (S == 3 | S == 4 | S == 6){
      fit.surv12 = cox.aalen(Surv(Ttilde1[delta1 == 1], Ttilde2[delta1 == 1], I(delta1*delta2)[delta1 == 1]) ~ 
                               #prop(L[delta1 == 1]) + prop(A[delta1 == 1]) + prop(I(L[delta1 == 1]^2)), 
                               prop(A[delta1 == 1]),
                             max.clust = n, data = df) 
      k12 = dim(fit.surv12$cum)[1]  
      
      g12.a0 = exp( fit.surv12$gamma[1] * rep(0, n) ) #+ fit.surv12$gamma[1] * L + fit.surv12$gamma[3] * I(L^2))
      g12.a1 = exp( fit.surv12$gamma[1] * rep(1, n) ) #+ fit.surv12$gamma[1] * L + fit.surv12$gamma[3] * I(L^2))
    } else {
      fit.surv12 = cox.aalen(Surv(Ttilde1[delta1 == 1], Ttilde2[delta1 == 1], I(delta1*delta2)[delta1 == 1]) ~ 
                               prop(L[delta1 == 1]) + prop(A[delta1 == 1]),
                             max.clust = n, data = df) 
      k12 = dim(fit.surv12$cum)[1]        
      
      g12.a0 = exp( fit.surv12$gamma[2] * 0 + fit.surv12$gamma[1] * L )
      g12.a1 = exp( fit.surv12$gamma[2] * 1 + fit.surv12$gamma[1] * L )
    }
    
    
    
    
    if (S == 2 | S == 4 | S == 6){
      fitA = glm(A ~ 1, binomial) # L + I(L^2)
    } else {
      fitA = glm(A ~ L, binomial)
    }
    pA1 = as.numeric(predict(fitA, type="response"))
    
    
    #### CENSORING COX ####
    
    # Just one model for censoring
    if (S == 5 | S == 6){
      
      fit.cens = cox.aalen(Surv(Ttilde2, (1-delta2)) ~ prop(A), max.clust = n, data = df) 
      kc = dim(fit.cens$cum)[1]  
      
      gc.a0 = exp( fit.cens$gamma[1] * rep(0, n))
      gc.a1 = exp( fit.cens$gamma[1] * rep(1, n))
      
    } else {
      
      fit.cens = cox.aalen(Surv(Ttilde2, (1-delta2)) ~ prop(A) + prop(1 * I(L > 0.5)), max.clust = n, data = df) 
      kc = dim(fit.cens$cum)[1]  
      
      gc.a0 = exp( fit.cens$gamma[1] * 0 + fit.cens$gamma[2] * 1 * I(L > 0.5))
      gc.a1 = exp( fit.cens$gamma[1] * 1 + fit.cens$gamma[2] * 1 * I(L > 0.5))
    }
    
    #### Various operations of changing support ####
    
    tau1 = fit.surv1$cum[,1]        # extracting the time vector (first column of fit.surv$cum)
    tau2 = fit.surv2$cum[,1]      
    tau12 = fit.surv12$cum[,1]   
    tauc = fit.cens$cum[,1] 
    
    #### Express everything wrt \Lambda_1 and compute d\Lambda_1 ####
    
    Lam1.0 = fit.surv1$cum[,2]
    Lam2.0.tau1 = apply_look_up_f_cpp(tau1, fit.surv2$cum[,1], fit.surv2$cum[,2])
    Lam12.0.tau1 = apply_look_up_f_cpp(tau1, fit.surv12$cum[,1], fit.surv12$cum[,2])
    
    Lam1.a0.tau1 = vvmult(g1.a0, Lam1.0)
    Lam2.a0.tau1 = vvmult(g2.a0, Lam2.0.tau1)
    Lam12.a0.tau1 = vvmult(g12.a0, Lam12.0.tau1)
    
    Lam1.a1.tau1 = vvmult(g1.a1, Lam1.0)
    Lam2.a1.tau1 = vvmult(g2.a1, Lam2.0.tau1)
    Lam12.a1.tau1 = vvmult(g12.a1, Lam12.0.tau1)
    
    dLam1.a0 = rbind(rep(0,n), Lam1.a0.tau1[2:k1,] - Lam1.a0.tau1[1:(k1-1),])
    dLam1.a1 = rbind(rep(0,n), Lam1.a1.tau1[2:k1,] - Lam1.a1.tau1[1:(k1-1),])
    
    
    #### Express everything wrt \Lambda_{12} and compute d\Lambda_{12} ####
    
    Lam12.0 = fit.surv12$cum[,2]
    Lam1.0.tau12 = apply_look_up_f_cpp(tau12, fit.surv1$cum[,1], fit.surv1$cum[,2])
    
    Lam12.a0.tau12 = vvmult(g12.a0, Lam12.0)
    Lam1.a0.tau12 = vvmult(g1.a0, Lam1.0.tau12)
    
    Lam12.a1.tau12 = vvmult(g12.a1, Lam12.0)
    Lam1.a1.tau12 = vvmult(g1.a1, Lam1.0.tau12)
    
    dLam12.a0 = rbind(rep(0,n), Lam12.a0.tau12[2:k12,] - Lam12.a0.tau12[1:(k12-1),])
    dLam12.a1 = rbind(rep(0,n), Lam12.a1.tau12[2:k12,] - Lam12.a1.tau12[1:(k12-1),])
    
    
    #### Express everything wrt \Lambda_{C} and compute d\Lambda_{C} ####
    
    LamC = fit.cens$cum[,2]
    Lam1.tauc = apply_look_up_f_cpp(tauc, fit.surv1$cum[,1], fit.surv1$cum[,2])
    Lam2.tauc = apply_look_up_f_cpp(tauc, fit.surv2$cum[,1], fit.surv2$cum[,2])
    Lam12.tauc = apply_look_up_f_cpp(tauc, fit.surv12$cum[,1], fit.surv12$cum[,2])
    
    LamC.a0.tauc = vvmult(gc.a0, LamC)
    LamC.a1.tauc = vvmult(gc.a1, LamC)
    
    Lam1.a0.tauc = vvmult(g1.a0, Lam1.tauc)
    Lam2.a0.tauc = vvmult(g2.a0, Lam2.tauc)
    Lam12.a0.tauc = vvmult(g12.a0, Lam12.tauc)
    
    Lam1.a1.tauc = vvmult(g1.a1, Lam1.tauc)
    Lam2.a1.tauc = vvmult(g2.a1, Lam2.tauc)
    Lam12.a1.tauc = vvmult(g12.a1, Lam12.tauc)
    
    KC.a0.tauc = exp(- LamC.a0.tauc )
    KC.a1.tauc = exp(- LamC.a1.tauc )
    
    dLamC.a0 = rbind(rep(0,n), LamC.a0.tauc[2:kc,] - LamC.a0.tauc[1:(kc-1),])
    dLamC.a1 = rbind(rep(0,n), LamC.a1.tauc[2:kc,] - LamC.a1.tauc[1:(kc-1),])
    
    
    #### One-step estimator ####
    
    estOS = OneStep(t, k, df, pA1,
                    tau12, dLam12.a0, dLam12.a1,
                    tau1, dLam1.a0, dLam1.a1,
                    Lam1.a0.tau1, Lam1.a1.tau1, 
                    Lam2.a0.tau1, Lam2.a1.tau1,
                    Lam1.a0.tauc, Lam1.a1.tauc, 
                    Lam2.a0.tauc, Lam2.a1.tauc,
                    Lam12.a0.tau1, Lam12.a1.tau1,
                    Lam12.a0.tau12, Lam12.a1.tau12, 
                    Lam12.a0.tauc, Lam12.a1.tauc,
                    tauc, KC.a0.tauc, KC.a1.tauc, dLamC.a0, dLamC.a1)
    
    return(estOS)
  })
}


#### One run ####

onerun = function(i, n = 1000, t = 10, k, S, CENS){
  set.seed(i)
  df = sim.data(sample.size = n, CENS)
  res = est.f(df, t, k, S) 
  return( res ) 
}



#### Compute coverage ####

coverage_function <- function(values, se_values, target_val, confidence_level = 0.95) {
  half_width <- qnorm(1 - (1 - confidence_level) / 2) * se_values
  lower <- values - half_width
  upper <- values + half_width
  coverage <- 1 * ((lower <= target_val) & (target_val <= upper))
  return(coverage)
}



#### Compute true value ####


truevalue = function(A, t = 10, k = 3, maxiter = 100000){
  set.seed(1)
  
  n = 1
  i = 1
  
  while( i <= maxiter){
    
    L = runif(n, min=0, max=1) # covariate
    #L = 1*I(L > 0.5)
    
    lambda1 = 0.04 * exp(log(2) * L + 1 * A)
    lambda2 = 0.02 * exp(log(2) * L + 1 * A)
    
    lambda12 = 0.05 * exp(log(2) * L + 1 * A)
    
    f <- function(t2, t1) ((exp(-lambda12*(t2-t1))) / pmin(t2, t)^(1/k)) * (exp(-t1*(lambda1 + lambda2))) * lambda12 * lambda1
    res = c(res, 
            integrate(Vectorize(function(t1) { 
              integrate(function(t2) { 
                f(t2, t1) 
              }, t1, Inf)$value
            }), 0, t)$value
    )
    i = i + 1
    set.seed(i)
  }
  return(mean(res))
}


#### Analyse results ####


ana <- function(res, TV1, TV0) {
  coef <- res[,1:2]
  scoef <- res[,3:4]
  plugin <- res[,5:6]
  debias <- res[,7:8]
  cov = numeric(ncol(coef))
  bias = numeric(2)
  for (i in seq(1, ncol(coef), by = 2)) {
    # Apply function to even columns
    cov[i] <- mean(coverage_function(coef[, i], scoef[, i], TV1 )) 
    bias[1] <- mean(coef[, i] - TV1) 
    # Apply function to odd columns
    cov[i+1] <- mean(coverage_function(coef[, i+1], scoef[, i+1], TV0 )) 
    bias[2] <- mean(coef[, i+1] - TV0) 
  }
  m <- round(cbind( apply(coef,2,mean), apply(plugin,2,mean), apply(debias,2,mean), bias, apply(coef,2,sd), apply(scoef,2,mean), cov), 3)
  colnames(m) <- c('One-step', 'Plug-in', 'De-bias', 'Bias', 'SD', 'SE', 'Cov')
  
  rownames(m) <- c("A = 1", "A = 0")
  
  return(round(m, 3))
}

