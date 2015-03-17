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
load("AllCov/pvalue05/Model1.Line.Diet.Blockorder.Block.RFI.RINb.Concb.RINa.Conca.neut.lymp.mono.eosi.baso/Model1_result.RData")
load("AllCov/pvalue05/Model1.Line.Diet.Blockorder.Block.RFI.RINb.Concb.RINa.Conca.neut.lymp.mono.eosi.baso/Model1_fit.RData")
dl2 <- model.matrix(~Line + Diet + Blockorder + Block + 
                      RFI + RINb + Concb + RINa + Conca + 
                      neut + lymp + mono + eosi + baso)

var.EBAllCov <- laply(1:dim(counts)[1], function(i){
  var.b1(b= fit$coef[i,], nbd = fit$NB.disp[i], dl2 = dl2)
})
# save(var.EBAllCov, file = "var.EBAllCov.RData")
load(file = "var.EBAllCov.RData")


ebhyperpar.EBAllCov <- laply(1:dim(dl2)[2], function(i) 
  EBhyperpar(var.EBAllCov[,i], fit$coef[,i]))

save(ebhyperpar.EBAllCov, file = "ebhyperpar.EBAllCov.RData")

 
coefEBAllCov <- laply(1:dim(fit$coef)[2], 
                      function(i) post.mu(ebhyperpar.EBAllCov[i,], var.EBAllCov[,i],
                                          fit$coef[,i]))
save(coefEBAllCov, file = "coefEBAllCov.RData")
# 
# 
# betavec <-t(dl2[,-c(1,2)]%*%coefEBAllCov [-c(1,2),])
# design.list.beta <- vector("list", 2)
# design.list.beta[[1]] <- model.matrix(~Line) # line
# design.list.beta[[2]] <-rep(1, 31)
# den.df <- dim(dl2)[2]
# fitEBAllCov <- myQL.fit(counts, design.list = design.list.beta, betavec = betavec,   # dim(counts)
#                                      log.offset = log.offset, 
#                                      print.progress=TRUE,
#                                      Model = "NegBin", method = "optim", den.df = den.df)
# resultEBAllCov <- QL.results(fitEBAllCov, Plot = F)
# save(fitEBAllCov, file = "fitEBAllCov.RData")
# save(resultEBAllCov, file = "resultEBAllCov.RData")