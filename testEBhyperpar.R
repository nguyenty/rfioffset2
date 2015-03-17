source("lf.R")
source("EBhyperpar.R")
library(plyr)
mu0 <- 0      
s0 <- .000001
mu1 <- 4
s1 <- 2
n <- 1000
pi0 <- .5
z <- rbinom(n, 1, pi0)
mu <- z*rnorm(n,mu0,s0) + (1-z)*rnorm(n,mu1,sqrt(s1))
sh <- runif(n, 0.04, 0.05)
y <- rnorm(n, mu, sqrt(sh))

EBhyperpar(sh, y)
