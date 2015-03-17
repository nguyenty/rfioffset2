
f <- function(m) class(try(solve(m),silent=T))=="matrix"

var.b1 <- function(b,nbd,dl2){
  # b is fitted coefficient, w is the negative binomial dispersion
  # dl2 is the design  matrix full model 
  # counts is the counts data for 1 gene
  X <- dl2
  
  eta <- c(log.offset + X%*%b)
  
  mu <- exp(eta)
  
  fish <- t(X)%*%(diag(nrow(X))*mu/(mu*nbd+1))%*%X
  
  if(f(fish)){
    fish.inv <- solve(fish)
    return(diag(fish.inv))
  }else{
    return(rep(0, ncol(X)))
  }
}
