source("simulations/utils.R")

library(xtable)

load_and_summarize1 <- function(varz, cens, which = c("Tab1","TabE1")) {
  
  which <- match.arg(which)
  
  file <- paste0(
    "simulations/intermediate_results/Tab1TabE1/", which,
    "_varz", varz,
    "_cens", cens,
    ".rda"
  )
  
  load(file)  # loads outtab1 or outtabE1 into environment
  
  OUTPUT <- if (which == "Tab1") outtab1 else outtabE1
  
  mm <- ana(OUTPUT$resl, true = OUTPUT$true)
  
  return(list(cbind(dep, varz, cens, mm$summary), cbind(dep, varz, cens, mm$coef)))
}


depvals <- c(1) # dependence: 1 = shared frailty (v=1), 4 = frailty only recurrent events (v=0)
typevals <- c(2)
cens_vals <- c(1/4,2/4)
theta_vals <- c(0.5,1,2)
scale1_vals <- c(1)
scaled_vals <- c(1)

dep = depvals[1]

for (which in c("Tab1", "TabE1")) {
  
  outtotSUMM <- outtotCOEF <- NULL
  
  for (varz in theta_vals) 
    for (cens in cens_vals) {
      
      row <- load_and_summarize1(
        varz   = varz,
        cens   = cens,
        which  = which
      )
      
      
      outtotSUMM <- rbind(outtotSUMM, row[[1]])
      outtotCOEF <- rbind(outtotCOEF, row[[2]])
    }
  
  outtotSUMM <- as.data.frame(outtotSUMM)
  outtotCOEF <- as.data.frame(outtotCOEF)
  
  ## Reproduce original column ordering
  temp <- NULL
  for (varz in theta_vals)
    for (cens in cens_vals) {
      df_row1 = outtotCOEF[(outtotCOEF$varz == varz) & 
                             (outtotCOEF$cens == cens), ][c(1,2,5,6),]
      rownames(df_row1) = NULL
      
      df_row <- outtotSUMM[(outtotSUMM$varz == varz) &
                             (outtotSUMM$cens == cens), -c(7)] # remove power from here
      rownames(df_row) = NULL
      
      cols = c(4,8,5,6,9,7)
      if (which == "Tab1"){
        temp = rbind(temp,
                     cbind(df_row1[1:2,c(2,3,cols)], df_row1[3:4,c(cols)]),
                     cbind(df_row[2,c(2,3,cols)], df_row[1,c(cols)]))
      } else {
        temp = rbind(temp,
                     df_row1[1:2,c(2,3,cols)],
                     df_row[2,c(2,3,cols)])
      }
      
      rownames(temp) = NULL
    }
  
  
  name <- if (which == "Tab1") "1" else "E1"
  
  print(xtable(temp, digits = 3))
  
  print(
    xtable(temp, digits = 3),
    file = paste0("simulations/results/Table", name, "Latex.txt"),
    include.rownames = FALSE
  )
}

