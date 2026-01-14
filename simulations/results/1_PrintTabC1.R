

# import functions
source("simulations/utilsTabC1.R")

# import results
load_res <- function(scen, cens) {
  
  file <- paste0(
    "simulations/intermediate_results/TabC1TabC2/C1results",
    "_", scen,
    "_", cens,
    "_insidek3.txt"
  )
  
  return(read.table(file, header=T))
}



# import results

# TRUE VALUES for cube root and t = 10
TV1 = 0.3423482  # truevalue(A = 1, t = 10, k = 3, maxiter = 1000000) 
TV0 = 0.1913787  # truevalue(A = 0, t = 10, k = 3, maxiter = 1000000) 




scenarios <- c(
  "(i) All correct",
  "(ii) Propen. score (PS) missp.",
  "(iii) $\\Lambda_{01}$ and $\\Lambda_{1D}$ missp.",
  "(iv) $\\Lambda_{01}$, $\\Lambda_{1D}$ and PS missp.",
  "(v) $\\Lambda_{\\widetilde{C}}$ missp.",
  "(vi) $\\Lambda_{\\widetilde{C}}$, $\\Lambda_{01}$, $\\Lambda_{1D}$ and PS missp."
)


alphas <- c(0.01, 0.03, 0.05)


results <- list()
k <- 1

for (s in 1:length(scenarios)) {
  results[[s]] <- list()
  
  for (a in alphas) {
    m <- ana(load_res(scen=s, cens=a), TV1, TV0)
    results[[s]][[paste0("alpha=", a)]] <- m
    k <- k + 1
  }
}


latex_rows <- function(m) {
  apply(m, 1, function(x) {
    paste(x, collapse = " & ")
  })
}


out_file <- "simulations/results/TableC1body.txt"
sink(out_file)

for (s in 1:length(scenarios)) {
  cat("\\multirow{6}{*}{", scenarios[s], "} \n", sep = "")
  
  for (a in names(results[[s]])) {
    cat("& \\multirow{2}{*}{ $\\alpha =", sub("alpha=", "", a), "$ } \n")
    
    rows <- latex_rows(results[[s]][[a]])
    cat("   & A = 1 &", rows[1], "\\\\\n")
    cat("&  & A = 0 &", rows[2], "\\\\\n")
    alpha_val <- sub("alpha=", "", a)
    if (alpha_val == "0.05") {
      cat("\\midrule\n")
    } else {
      cat("\\cmidrule{2-10}\n")
    }
  }
}


