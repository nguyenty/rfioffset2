source("lf.R")
EBhyperpar <- function(sh, coef){
#   par0 <- c(.5, mean(coef), var(coef))
  par0 <- c(.8, mean(coef)*5, var(coef))
  outpar <- optim(par = par0, lf, method = "BFGS", 
                  lower = c(0, -Inf, 0), upper = c(1, Inf, 100), sh = sh, 
                  coef = coef)
  par <- outpar$par
  par
}