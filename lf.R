
lf <- function(par, sh,  coef){
  ys <- coef
  l <-  laply(1:length(ys), function(i)
    par[1]*dnorm(ys[i], 0, sqrt(.000001^2+ sh[i])) + 
      (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + sh[i])))
  -sum(log(l))
}