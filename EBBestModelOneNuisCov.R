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

load("BestModelOneNuisCov/pvalue05/Model1.Line.Block.Concb.RINa.neut.lymp.mono.eosi.baso/Model1_result.RData")
load("BestModelOneNuisCov/pvalue05/Model1.Line.Block.Concb.RINa.neut.lymp.mono.eosi.baso/Model1_fit.RData")
dl2 <- model.matrix(~Line + Block + Concb + RINa + neut + lymp + mono + eosi + baso)

var.EBBestModelOneNuisCov <- laply(1:dim(counts)[1], function(i){
  var.b1(b= fit$coef[i,], w = fit$NB.disp[i], dl2 = dl2)
})
save(var.EBBestModelOneNuisCov, file = "var.EBBestModelOneNuisCov.RData")
load(file = "var.EBBestModelOneNuisCov.RData")

ebhyperpar.EBBestModelOneNuisCov <- laply(1:dim(fit$coef)[2], 
                                          function(i) 
                                            EBhyperpar(i, var.EBBestModelOneNuisCov))
save(ebhyperpar.EBBestModelOneNuisCov, file = "ebhyperpar.EBBestModelOneNuisCov.RData")
load(file = "ebhyperpar.EBBestModelOneNuisCov.RData")

coefEBBestModelOneNuisCov <- laply(1:dim(fit$coef)[2], function(i) 
  post.mu(i,ebhyperpar.EBBestModelOneNuisCov, var.EBBestModelOneNuisCov))
save(coefEBBestModelOneNuisCov, file = "coefEBBestModelOneNuisCov.RData")



load(file = "coefEBBestModelOneNuisCov.RData")

betavec <-t(dl2[,-c(1,2)]%*%coefEBBestModelOneNuisCov [-c(1,2),])

design.list.beta <- vector("list", 2)
design.list.beta[[1]] <- model.matrix(~Line) # line
design.list.beta[[2]] <-rep(1, 31)
den.df <- dim(dl2)[2]
fitEBBestModelOneNuisCov <- myQL.fit(counts, design.list = design.list.beta, betavec = betavec,   # dim(counts)
                    log.offset = log.offset, 
                    print.progress=TRUE,
                    Model = "NegBin", method = "optim", den.df = den.df)


resultEBBestModelOneNuisCov <- QL.results(fitEBBestModelOneNuisCov, Plot = F)

hist(resultEBBestModelOneNuisCov$P.values[[3]])
save(fitEBBestModelOneNuisCov, file = "fitEBBestModelOneNuisCov.RData")
save(resultEBBestModelOneNuisCov, file = "resultEBBestModelOneNuisCov.RData")
