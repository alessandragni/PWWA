#####  Case Study based on HF-Action trial - Section 6 #####
# Code used for producing Figures 2, 3 and Table 3 in Section 6

###### Import packages ######
library(mets) 
library(latex2exp)

###### Load data and preprocessing ######

load("data/hhosp.rda")
dfactor(hhosp) <- treatment~trt_ab
hhosp$status13 <- hhosp$status
hhosp <- dtransform(hhosp,status13=1,status==3)
aahosp <- dkeep(hhosp,~start_time+stop_time+status+status13+etiology+priormi+treatment+trt_ab+sex+region+age+beckb+hfhosp+id)
rr <- na.omit(aahosp)
rr <- count.history(rr)
rrid <- countID(rr)
rr <- cbind(rr,rrid)
rr$time0 <- rr$start_time
rr$time1 <- rr$stop_time
## only look at etiology==2
rr <- subset(rr,etiology==2)
rr <- cbind(countID(rr),rr)
## let id be 1, ..., as data appears
rr$id <- rr$indexid+1



###### Data for Figures 2 and 3 ######

# Cycle to fill data of to create plots in Section 6
res <- resR <- resNUM <- resDEN <- c()
for (tt in seq(0.5,3.9,by=0.1)) {
  dd <- WA_recurrent(Event(time0,time1,status)~treatment+cluster(id),rr,time=tt,death.code=3,
                      augmentR=~sex+region+age+beckb+hfhosp,augmentC=~sex+region+age+beckb+hfhosp)
  
  # PWWA
  pwwa <- estimate(coef=dd$ET$riskDR$riskDR,vcov=dd$ET$riskDR$var.riskDR) 
  pwwa_diff = estimate(pwwa,contrast=rbind(c(1,-1)))
  
  # EWWA
  ratio <- estimate(dd$RAW$ratio.means) 
  ratio_diff <- estimate(ratio,contrast=rbind(c(1,-1)))
  
  # NUM
  num <- estimate(dd$RAW$meanN)
  num_diff <- estimate(num,contrast=rbind(c(1,-1)))
  
  # DENOM
  den <- estimate(dd$RAW$rmst)
  den_diff <- estimate(den,contrast=rbind(c(1,-1)))
  
  # rbind of time, estimates (0 and 1), lower CI, upper CI, p-value of difference
  res <- rbind(res,c(tt,pwwa$coef,pwwa$coefmat[1,3:4],pwwa$coefmat[2,3:4],pwwa_diff$coefmat[,5]))
  resR <- rbind(resR,c(tt,ratio$coef,ratio$coefmat[1,3:4],ratio$coefmat[2,3:4],ratio_diff$coefmat[,5]))
  resNUM <- rbind(resNUM, c(tt, num$coef, num$coefmat[1,3:4], num$coefmat[2,3:4], num_diff$coefmat[,5]))
  resDEN <- rbind(resDEN, c(tt, den$coef, den$coefmat[1,3:4], den$coefmat[2,3:4], den_diff$coefmat[,5]))
}


# Function to create plots of Section 6
plotres <- function(res, num = FALSE, ylab) {
  if(num==TRUE){
    matplot(res[,1],res[,2:3],type="l",lwd=3,ylim=c(0,3),xlab="Time (years)",ylab=ylab)
  } else {
    matplot(res[,1],res[,2:3],type="l",lwd=3,ylim=c(0.5,1.2),xlab="Time (years)",ylab=ylab)
  }
  plotConfRegion(res[,1],res[,4:5],col=1)
  plotConfRegion(res[,1],res[,6:7],col=2)
  if (ncol(res)==8) {
    sigp <- (res[,8]<0.05)
    points(res[sigp,1],res[sigp,2],pch="*",cex = 2.5, col=1)
    points(res[sigp,1],res[sigp,3],pch="*",cex = 2.5, col=2)
  }
  legend("bottomright",c("usual care","training"),lty=1:2,col=1:2,lwd=2.5)
}


###### Figure 2 - Section 6 ######
pdf("casestudies/results/Fig2-hf-action-wwa-NUM.pdf", width=6, height=5)
plotres(resNUM, num = TRUE, ylab = TeX("$E[N(t \\wedge T_D)]$"))
dev.off()

###### Figure 3 (i) - Section 6 ######
pdf("casestudies/results/Fig3i-hf-action-wwa-PWWA.pdf", width=5, height=5)
plotres(res, ylab = "PWWA estimand")
dev.off()

###### Figure 3 (ii) - Section 6 ######
pdf("casestudies/results/Fig3ii-hf-action-wwa-EWWA.pdf", width=5, height=5)
plotres(resR, ylab = "EWWA estimand")
dev.off()







###### Table 3 - Section 6 ######

#  vector of values
tt_values <- c(2.3, 3.3, 3.9)

# Run WA_recurrent and store results
results_list <- list()
for(tt in tt_values){ 
  results_list[[as.character(tt)]] <- WA_recurrent(
    Event(time0,time1,status)~treatment+cluster(id),
    rr, time=tt, death.code=3,
    augmentR=~sex+region+age+beckb+hfhosp,
    augmentC=~sex+region+age+beckb+hfhosp
  )
}

# Function to extract estimates, SE, and p-values
extract_table <- function(dd){

  num <- estimate(dd$RAW$meanN)
  denom <- estimate(dd$RAW$rmst)
  ratio <- estimate(dd$RAW$ratio.means)
  
  num_diff <- estimate(num,contrast=rbind(c(1,-1)))
  denom_diff <- estimate(denom,contrast=rbind(c(1,-1)))
  ratio_diff <- estimate(ratio,contrast=rbind(c(1,-1)))
  
  list(
    numerator = paste0(round(num$coefmat[,1],3), " (", round(num$coefmat[,2],3), ")"),
    numerator_p = round(num$coefmat[,5],3),
    
    numerator_diff = paste0(round(num_diff$coefmat[,1],3), " (", round(num_diff$coefmat[,2],3), ")"),
    numerator_diff_p = round(num_diff$coefmat[,5],3),
    
    denominator = paste0(round(denom$coefmat[,1],3), " (", round(denom$coefmat[,2],3), ")"), 
    denominator_p = round(denom$coefmat[,5],3),
    
    denominator_diff = paste0(round(denom_diff$coefmat[,1],3), " (", round(denom_diff$coefmat[,2],3), ")"), 
    denominator_diff_p = round(denom_diff$coefmat[,5],3),
    
    ratio_diff = paste0(round(ratio_diff$coefmat[,1],3), " (", round(ratio_diff$coefmat[,2],3), ")"),
    ratio_p = round(ratio_diff$coefmat[,5],3)
  )
}


# Function to generate latex table
generate_latex_table <- function(results_list, tt_values) {
  
  # --- Round p-values < 0.001 to 0
  format_p <- function(p) ifelse(p < 0.001, 0, signif(p, 3))
  
  # --- Initialize vectors for rows
  num_treat0 <- num_treat1 <- num_diff <- c()
  num_treat0_p <- num_treat1_p <- num_diff_p <- c()
  
  den_treat0 <- den_treat1 <- den_diff <- c()
  den_treat0_p <- den_treat1_p <- den_diff_p <- c()
  
  ratio_diff <- c()
  ratio_p <- c()
  
  # --- Fill vectors row by row
  for(tt in tt_values){
    dd <- extract_table(results_list[[as.character(tt)]])
    
    # Numerator
    num_treat0 <- c(num_treat0, dd$numerator[1])
    num_treat1 <- c(num_treat1, dd$numerator[2])
    num_diff   <- c(num_diff, dd$numerator_diff)
    
    num_treat0_p <- c(num_treat0_p, format_p(dd$numerator_p[1]))
    num_treat1_p <- c(num_treat1_p, format_p(dd$numerator_p[2]))
    num_diff_p   <- c(num_diff_p, format_p(dd$numerator_diff_p))
    
    # Denominator
    den_treat0 <- c(den_treat0, dd$denominator[1])
    den_treat1 <- c(den_treat1, dd$denominator[2])
    den_diff   <- c(den_diff, dd$denominator_diff)
    
    den_treat0_p <- c(den_treat0_p, format_p(dd$denominator_p[1]))
    den_treat1_p <- c(den_treat1_p, format_p(dd$denominator_p[2]))
    den_diff_p   <- c(den_diff_p, format_p(dd$denominator_diff_p))
    
    # Ratio
    ratio_diff <- c(ratio_diff, dd$ratio_diff)
    ratio_p    <- c(ratio_p, format_p(dd$ratio_p))
  }
  
  # --- Build a LaTeX row
  latex_row <- function(est, pval) {
    paste0(est[1], " & ", pval[1], " & ",
           est[2], " & ", pval[2], " & ",
           est[3], " & ", pval[3], " \\\\")
  }
  
  # --- Build LaTeX table as a string
  table_code <- paste0(
    "\\begin{table}[H]\n\\centering \n\\resizebox{\\textwidth}{!}{%\n",
    "\\begin{tabular}{@{}llcccccc@{}}\n\\toprule\n",
    "& & \\multicolumn{2}{c}{$t=", tt_values[1], "$} & \\multicolumn{2}{c}{$t=", tt_values[2], "$} & \\multicolumn{2}{c}{$t=", tt_values[3], "$} \\\\\n",
    "\\cmidrule(l{3pt}r{3pt}){3-4} \\cmidrule(l{3pt}r{3pt}){5-6} \\cmidrule(l{3pt}r{3pt}){7-8}\n",
    "& & Est (SE) & $p$-val & Est (SE) & $p$-val & Est (SE) & $p$-val \\\\\n\\midrule\n",
    "\\multirow{3}{*}{Numerator}\n",
    "& usual care (0) & ", latex_row(num_treat0, num_treat0_p), "\n",
    "& training (1) & ", latex_row(num_treat1, num_treat1_p), "\n\\cmidrule{2-8}\n",
    "& Difference (0 - 1) & ", latex_row(num_diff, num_diff_p), "\n\\midrule\n",
    "\\multirow{3}{*}{Denominator}\n",
    "& usual care (0) & ", latex_row(den_treat0, den_treat0_p), "\n",
    "& training (1) & ", latex_row(den_treat1, den_treat1_p), "\n\\cmidrule{2-8}\n",
    "& Difference (0 - 1) & ", latex_row(den_diff, den_diff_p), "\n\\midrule\n",
    "Ratio (EWWA) & Difference (0 - 1) & ", latex_row(ratio_diff, ratio_p), "\n",
    "\\bottomrule\n\\end{tabular}}\n",
    "\\end{table}"
  )
  
  return(table_code)
}


# Generate LaTeX table code for Table 3
latex_code <- generate_latex_table(results_list, tt_values)
cat(latex_code)


# Define output path
out_file <- "casestudies/results/Table3Latex.txt"

# Write LaTeX code to txt file
writeLines(latex_code, con = out_file)

