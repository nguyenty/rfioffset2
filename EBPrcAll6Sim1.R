library("plyr")
source("readdata.R")
source("design_list.R")
source("fit_model.R")
source("simulation_criteria.R")
source("quasiseq shrinkage functions3.R")
source("varCoef.R")# b,nbd,dl2
source("lf.R") # par, sh,  coef
source("EBhyperpar.R") # sh, coef
source("post.mu.R") # hyperpar, sh, coef

#load("~/research/rfi/pvalue0511/pi0_0.6/simdat_1.RData")
pm1 <- proc.time()
vnf <- c("Line", "Diet", "Blockorder", "Block", "RFI", "RINb", "Concb", 
         "RINa", "Conca", "neut", "lymp", 
         "mono", "eosi", "baso")
fm <- as.formula(paste(" ~  ", paste(vnf[-1], collapse = "+")))
full_model <- scale(model.matrix(fm)[,-1])
u <- svd(full_model)$u
# dim(u)
colnames(u) <- paste("Prc", 1:dim(u)[2], sep = "")
for (i in 1:dim(u)[2]) assign(paste("Prc", i, sep =''), u[,i])

dl2 <- model.matrix(as.formula(
  paste(" ~ ", paste(vnf[1], collapse = "+"), 
        "+", paste(colnames(u), collapse = "+"))
))

vnf2 <- c(vnf[1], colnames(u))
vntest2 <- vnf[1]
list_out <- list_model(vnf2, vntest2)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
#dir <- "/run/user/1000/gvfs/smb-share:server=smb.stat.iastate.edu,share=ntyet"
# dir <- "~"
R6sim1 <- S6sim1 <- FDR6sim1 <- ind6sim1 <- NULL
result6sim1 <- list()
# resultsim1[[1]] <- 1
# resultsim1[[2]] <- 2
for (i in 1:100){# i <- 2
#   print(paste("i= ", i))
  path <- paste0("~/research/rfi/pvalue0511/pi0_0.6/simdat_", i,".RData")
  #  path <- paste0("P:/research/rfi/pvalue0511/pi0_0.6/simdat_", i,".RData")
  #  path <- paste0("/Volumes/ntyet/research/rfi/pvalue0511/pi0_0.6/simdat_", i,".RData")
  load(path)
  y <- simdat$y
  log.offset <- log(apply(y, 2, quantile, 0.75))
  
  
  fit <- QL.fit(y, design.list, test.mat, # dim(counts)
                log.offset = log.offset, print.progress=FALSE,
                Model = "NegBin", method = "optim")
  result<- QL.results(fit, Plot = FALSE)
  
  var.EBPrcAllsim1 <- laply(1:dim(y)[1], function(j){
    var.b1(b= fit$coef[j,], nbd = fit$NB.disp[j], dl2 = dl2)
  })
  
  if (any(var.EBPrcAllsim1 ==0)){
    ind6sim1[i] <- 0
    print(paste("i = ", i, ", ind6sim1 = ", ind6sim1[i]))
  } else{
    ind6sim1[i] <- 1
    print(paste("i = ", i, ", ind6sim1 = ", ind6sim1[i]))
    ebhyperpar.EBPrcAllsim1 <- laply(1:dim(fit$coef)[2], function(j) 
      EBhyperpar(var.EBPrcAllsim1[,j], fit$coef[,j]))
    
#     ebhyperpar.EBPrcAllsim1 <- laply(1:11, function(j) 
#       EBhyperpar(var.EBPrcAllsim1[,j], fit$coef[,j]))
#     
    # save(ebhyperpar.EBPrcAllsim1, file = "ebhyperpar.EBPrcAllsim1.RData")
    
    coefEBPrcAllsim1 <- laply(1:dim(fit$coef)[2], function(j) 
      post.mu(ebhyperpar.EBPrcAllsim1[j,], var.EBPrcAllsim1[,j], fit$coef[,j]))
    # save(coefEBPrcAllsim1, file = "coefEBPrcAllsim1.RData")
    
    
    betavec <-t(dl2[,-c(1,2)]%*%coefEBPrcAllsim1 [-c(1,2),])
    design.list.beta <- vector("list", 2)
    design.list.beta[[1]] <- model.matrix(~Line) # line
    design.list.beta[[2]] <-rep(1, 31)
    den.df <-31-dim(dl2)[2]
    fitEBPrcAllsim1 <- myQL.fit(y, design.list = design.list.beta, betavec = betavec,   # dim(counts)
                                log.offset = log.offset, 
                                print.progress=FALSE,
                                Model = "NegBin", method = "optim", den.df = den.df)
    resultEBPrcAllsim1 <- QL.results(fitEBPrcAllsim1, Plot = F)
    # save(fitEBPrcAllsim1, file = "fitEBPrcAllsim1.RData")
    # save(resultEBPrcAllsim1, file = "resultEBPrcAllsim1.RData")
    #hist(resultEBPrcAllsim1$P.values[[3]], nclass = 100)
    R6sim1[i] <- sum(resultEBPrcAllsim1$Q.values[[3]] <=.05)
    S6sim1[i] <- sum(which(resultEBPrcAllsim1$Q.values[[3]] <=.05)>simdat$EE)
    FDR6sim1[i] <- 1 - S6sim1[i]/max(R6sim1[i],1)
    print(paste0('FDR= ', FDR6sim1[i], ", S = ", S6sim1[i], ", R = ", R6sim1[i]))
    result6sim1[[i]] <- resultEBPrcAllsim1
  }
  
}
save(result6sim1, file = "result6sim1.RData")
save(R6sim1, file = "R6sim1.RData")
save(S6sim1, file = "S6sim1.RData")
save(FDR6sim1, file = "FDR6sim1.RData")
save(ind6sim1, file = "ind6sim1.RData")
proc.time() -pm1