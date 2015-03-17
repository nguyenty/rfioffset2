library(fdrtool)

kss <- function(z){ # z <- result$P.values[[3]][,"Line"]
  e <- ecdf(z)
  g <- grenander(e)
  if(tail(g$x.knots, 1)!=1){
    gx <- c(0, g$x.knots, 1); gF <- c(0, g$F.knots, 1)
  }
  if (tail(g$x.knots, 1)==1){
    gx <- c(0, g$x.knots); gF <- c(0, g$F.knots)
  }
  out <- max(abs(gF - gx))
  return(out)  
} 

#load("U:/R/RA/Data/RFI-newdata/resultpairedcbc/ks/Model2.Line.Diet.RFI.Concb.Conca.RINa.neut.lymp.mono.eosi.baso.Block.Blockorder/Model2_result.RData")
# sel_criteria(result)
#load("U:/R/RA/Data/RFI-newdata/resultpairedcbcQLfit2Diet2/pvalue05/Model2.Line.Concb.RINb.Conca.RINa.neut.lymp.mono.baso.Block/Model2_result.RData")
sel_criteria <- function(result){ # sel_criteria(result)
  dat <- result$P.values[[3]][,colnames(result$P.values[[3]])]
  if(is.vector(dat)) dat <- as.matrix(dat)
   # Kolmogorow Smirnov statistics 
  ks <- apply(dat, 2, function(z)kss(z))
  pvalue_05 <- apply(dat<=0.05, 2, sum)
  out <- data.frame(pvalue05 = order(pvalue_05),
                    ks = order(ks))
  return(out)
}
#load("U:/R/RA/Data/RFI-newdata/resultpairedcbc/ks/Model8.Line.Concb.RINa.neut.lymp.mono.Block/Model8_result.RData")

#sel_criteria(result)
pauc_out <- function(test.vector, lab){
  lab <- as.factor(lab)
  roc.out <- roc(1-test.vector, lab) # plot(roc.out)
  roc.ind <- sum(roc.out$fpr<=.05)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- auc(roc.out, min =roc.min)
  return(pauc)
}
