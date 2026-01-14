

# import functions
source("simulations/utilsTabC2.R")

# import results
load_res <- function(scen = 1, cens) {
  
  file <- paste0(
    "simulations/intermediate_results/TabC1TabC2/C2results",
    "_", scen,
    "_", cens,
    "_insidek3.txt"
  )
  
  return(read.table(file, header=T))
}


alphas <- c(0.01, 0.03, 0.05)

datalist <- lapply(alphas, function(a) {
  load_res(cens = a)
})

# Generate LaTeX table body for each result
table_body <- mapply(function(data, alpha) {
  latex_body(data, alpha)
}, datalist, alphas, SIMPLIFY = TRUE)

# Print to console
cat(paste(table_body, collapse = "\n"))

# Write to file
writeLines(table_body, "simulations/results/TableC2body.txt")
