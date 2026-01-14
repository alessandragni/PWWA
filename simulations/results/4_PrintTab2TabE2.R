source("simulations/utils.R")

library(xtable)

load_and_summarize2 <- function(type, dep, varz, scaled, which = c("Tab2","TabE2")) {
  
  which <- match.arg(which)
  
  file <- paste0(
    "simulations/intermediate_results/Tab2TabE2/", which,
    "_type", type,
    "_dep", dep,
    "_varz", varz,
    "_scaled", scaled,
    ".rda"
  )
  
  load(file)  # loads outtab2 or outtabE2 into environment
  
  OUTPUT <- if (which == "Tab2") outtab2 else outtabE2
  
  mm <- anaR(OUTPUT$resl, true = OUTPUT$true)
  
  cbind(type, dep, varz, scaled, mm$summary)
}


depvals    <- c(1,4)
typevals   <- c(2,3,1)
theta_vals <- c(1,2)
scaled_vals<- c(1,4)


for (which in c("Tab2", "TabE2")) {
  
  outtotSUMM <- NULL
  
  for (type in typevals)
    for (dep in depvals)
      for (varz in theta_vals)
        for (scaled in scaled_vals) {
          
          row <- load_and_summarize2(
            type   = type,
            dep    = dep,
            varz   = varz,
            scaled = scaled,
            which  = which
          )
          
          outtotSUMM <- rbind(outtotSUMM, row)
        }
  
  outtotSUMM <- as.data.frame(outtotSUMM)
  
  ## Reproduce original column ordering
  temp <- NULL
  for (type in typevals)
    for (dep in depvals)
      for (varz in theta_vals)
        for (scaled in scaled_vals) {
          
          df_row <- outtotSUMM[
            outtotSUMM$type == type &
              outtotSUMM$dep == dep &
              outtotSUMM$varz == varz &
              outtotSUMM$scaled == scaled, ]
          
          if (nrow(df_row) == 0) next
          
          cols <- c(9,6,7,10,8) # 5 for mean
          temp <- rbind(
            temp,
            cbind(df_row[2, c(1:4, cols)], df_row[3, cols])
          )
        }
  
  name <- if (which == "Tab2") "2" else "E2"
  
  print(xtable(temp, digits = 3))
  
  print(
    xtable(temp, digits = 3),
    file = paste0("simulations/results/Table", name, "Latex.txt"),
    include.rownames = FALSE
  )
}
