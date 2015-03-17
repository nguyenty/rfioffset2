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

load("PrcAll/pvalue05/Model1.Line.Prc1.Prc2.Prc3.Prc4.Prc5.Prc6.Prc7.Prc8.Prc9.Prc10.Prc11.Prc12.Prc13.Prc14.Prc15.Prc16.Prc17.Prc18.Prc19.Prc20.Prc21/Model1_result.RData")
load("PrcAll/pvalue05/Model1.Line.Prc1.Prc2.Prc3.Prc4.Prc5.Prc6.Prc7.Prc8.Prc9.Prc10.Prc11.Prc12.Prc13.Prc14.Prc15.Prc16.Prc17.Prc18.Prc19.Prc20.Prc21/Model1_fit.RData")

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

var.EBPrcAll <- laply(1:dim(counts)[1], function(i){
  var.b1(b= fit$coef[i,], w = fit$NB.disp[i], dl2 = dl2)
})

save(var.EBPrcAll, file = "var.EBPrcAll.RData")
load(file = "var.EBPrcAll.RData")
ebhyperpar.EBPrcAll <- laply(1:dim(fit$coef)[2], function(i) 
  EBhyperpar(i, var.EBPrcAll, fit$coef))

save(ebhyperpar.EBPrcAll, file = "ebhyperpar.EBPrcAll.RData")

coefEBPrcAll <- laply(1:dim(fit$coef)[2], function(i) 
  post.mu(i,ebhyperpar.EBPrcAll, var.EBPrcAll, fit$coef))
save(coefEBPrcAll, file = "coefEBPrcAll.RData")


betavec <-t(dl2[,-c(1,2)]%*%coefEBPrcAll [-c(1,2),])
design.list.beta <- vector("list", 2)
design.list.beta[[1]] <- model.matrix(~Line) # line
design.list.beta[[2]] <-rep(1, 31)
den.df <- dim(dl2)[2]
fitEBPrcAll <- myQL.fit(counts, design.list = design.list.beta, betavec = betavec,   # dim(counts)
                                     log.offset = log.offset, 
                                     print.progress=TRUE,
                                     Model = "NegBin", method = "optim", den.df = den.df)
resultEBPrcAll <- QL.results(fitEBPrcAll, Plot = F)
save(fitEBPrcAll, file = "fitEBPrcAll.RData")
save(resultEBPrcAll, file = "resultEBPrcAll.RData")