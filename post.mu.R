post.mu <- function(hyperpar, sh, coef){
  ys <- coef
#   par <- c(.5, .1, 1)
#   outpar <- optim(par = par, lf, method = "BFGS", lower = c(0, -Inf, 0), upper = c(1, Inf, 100), index = index)
#   par <- outpar$par
  par <- hyperpar
  #if (par[3] == 0) par[3] <- .000001^2
  py <- laply(1:length(ys), function(i)
   par[1]*dnorm(ys[i], 0, sqrt(.000001^2+ sh[i])) + 
     (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + sh[i])))
   
   emu <- laply(1:length(ys), function(i)
     (par[1]*dnorm(ys[i], 0, sqrt(0.000001^2+sh[i]))* 
                (1/sh[i]*ys[i]+1/0.000001^2*0)/(1/sh[i]+1/0.000001^2)+ 
                (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + sh[i]))* 
                (1/sh[i]*ys[i]+1/par[3]*par[2])/(1/sh[i]+1/par[3]))/py[i])
  emu
  }
