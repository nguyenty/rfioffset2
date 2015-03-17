library("plyr")
source("readdata.R")
source("design_list.R")
source("fit_model.R")
source("simulation_criteria.R")
source("quasiseq shrinkage functions3.R")
source("varCoef.R")
source("lf.R")
source("EBhyperpar.R")
source("post.mu.R")
load(file = "var.EBPrcAll.RData")
n <- dim(var.EBPrcAll)[1]
pi <- .9
z <- rbinom(n, 1, pi)
mu0 <- 2
s20 <- 1
mui <- z*0 + (1-z)*rnorm(n, mu0, s20)
y <- rnorm(n, mui, var.EBPrcAll[,4])

pm <- proc.time()
EBhyperpar( var.EBPrcAll[,4],y )
proc.time() -pm
