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
load("Prc/pvalue05/Model1.Line.Diet.Blockorder.Block.Prc1.Prc2.Prc3.Prc4.Prc5.Prc6.Prc7.Prc8.Prc9.Prc10/Model1_result.RData")
load("Prc/pvalue05/Model1.Line.Diet.Blockorder.Block.Prc1.Prc2.Prc3.Prc4.Prc5.Prc6.Prc7.Prc8.Prc9.Prc10/Model1_fit.RData")
vnf <- c("Line", "Diet", "Blockorder", "Block", "RFI", "RINb", "Concb", 
         "RINa", "Conca", "neut", "lymp", 
         "mono", "eosi", "baso")
fm <- as.formula(paste(" ~ 0 + ", paste(vnf[-c(1:4)], collapse = "+")))
full_model <- scale(model.matrix(fm))
u <- svd(full_model)$u
colnames(u) <- paste("Prc", 1:10, sep = "")
for (i in 1:dim(u)[2]) assign(paste("Prc", i, sep =''), u[,i])

dl2 <- model.matrix(as.formula(
  paste(" ~ ", paste(vnf[c(1:4)], collapse = "+"), 
        "+", paste(colnames(u), collapse = "+"))
))


var.EBPrc <- laply(1:dim(counts)[1], function(i){
  var.b1(b= fit$coef[i,], w = fit$NB.disp[i], dl2 = dl2)
})

save(var.EBPrc, file = "var.EBPrc.RData")
load(file = "var.EBPrc.RData")
ebhyperpar.EBPrc <- laply(1:dim(fit$coef)[2], function(i) 
  EBhyperpar(i, var.EBPrc, fit$coef))

save(ebhyperpar.EBPrc, file = "ebhyperpar.EBPrc.RData")

coefEBPrc <- laply(1:dim(fit$coef)[2], function(i) 
  post.mu(i,ebhyperpar.EBPrc, var.EBPrc,fit$coef))
save(coefEBPrc, file = "coefEBPrc.RData")


betavec <-t(dl2[,-c(1,2)]%*%coefEBPrc [-c(1,2),])
design.list.beta <- vector("list", 2)
design.list.beta[[1]] <- model.matrix(~Line) # line
design.list.beta[[2]] <-rep(1, 31)
den.df <- dim(dl2)[2]
fitEBPrc <- myQL.fit(counts, design.list = design.list.beta, betavec = betavec,   # dim(counts)
                                     log.offset = log.offset, 
                                     print.progress=TRUE,
                                     Model = "NegBin", method = "optim", den.df = den.df)
resultEBPrc <- QL.results(fitEBPrc, Plot = F)
save(fitEBPrc, file = "fitEBPrc.RData")
save(resultEBPrc, file = "resultEBPrc.RData")