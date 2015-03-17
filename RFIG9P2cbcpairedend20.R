# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(Matrix)
library(edgeR)
library(reshape)
library(plyr)
library(fields)
library(reshape)
library(fdrtool)
source("QL.fit2.R")
source("QL.results.R")
source("NBDev4.R")
source("PoisDev.R")
source("negbin.br.R")
source("fbrglm.R")
source("glm.fit0.R")
library(AUC)

resultdir <- paste0(getwd(),"/resultpairedcbcQLfit20")
dir.create(resultdir, showWarnings = FALSE)
scount <- read.table("paired end uniquely mapped reads count table.txt", 
                     header = T)
scount <- scount[-c(which(scount[,1] %in%"ENSSSCG00000007978"),
                    which(scount[,1] %in%"ENSSSCG00000014725")),]
counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                             rowMeans(scount[,-1])>20 ,-1])

dim(counts)
###List of models function ####
covset <- read.table("covset.txt")
Blockorder <- as.factor(covset$Blockorder)
Block <- as.factor(covset$Block)
Line <- as.factor(covset$Line)
Diet <- as.factor(covset$Diet)
#table(Line, Diet)

RFI <- covset$RFI
RINa <- covset$RINa
RINb <- covset$RINb
Conca <- covset$Conca
Concb <- covset$Concb
neut <- covset$neut
lymp <- covset$lymp
mono <- covset$mono
baso <- covset$baso
eosi <- covset$eosi

#write.csv(counts, file = "counts_12280.csv")

log.offset <- log(apply(counts, 2, quantile, .75))

# str(fit)
# load("Model1_result.RData")
# sel_criteria(result)
#source("simulation_loadModel7.R")
fit_model <- function(full_model, model_th, criteria){ # model_th <- 1 criteria <- 1
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit <- QL.fit(counts, design.list, test.mat, # dim(counts)
                log.offset = log.offset, 
                print.progress=TRUE,
                Model = "NegBin", method = "optim")
  
  result<- QL.results(fit, Plot = FALSE)
  res_sel <- sel_criteria(result)
  k <- nrow(test.mat)
  name_model <- NULL 
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(resultdir,"/", colnames(res_sel)[criteria], "/Model",model_th,name_model, sep ="")
  dir.create(model_dir, showWarnings = FALSE)
  save(result, file = paste(model_dir,"/Model",model_th, "_result.RData", sep =""))
  save(fit, file = paste(model_dir,"/Model",model_th, "_fit.RData", sep =""))
  for(i in 1:(nrow(test.mat))){
    postscript(paste(model_dir,"/Model", 
                     model_th, row.names(test.mat)[i],".eps", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
    
    pdf(paste(model_dir,"/Model", 
              model_th, row.names(test.mat)[i],".pdf", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
  }
  print(paste("Model", model_th, sep = " "))
  
  return(res_sel)
}


source("simulation_criteria.R")
source("simulation_listmodel.R")
# load("Model1_result.RData")
# 
# load("U:/R/RA/Data/RFI-newdata/resultpairedcbc/ad/Model2.Line.Diet.RFI.Concb.RINb.RINa.neut.lymp.mono.eosi.baso.Block.Blockorder/Model2_result.RData")
# load("U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue05/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.neut.lymp.mono.eosi.baso.Block.Blockorder/Model1_result.RData")
# load("U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue05/Model2.Line.Diet.RFI.Concb.Conca.RINa.neut.lymp.mono.eosi.baso.Block.Blockorder/Model2_result.RData")
# load("U:/R/RA/Data/RFI-newdata/resultpairedcbc/ad/Model2.Line.Diet.RFI.Concb.RINb.RINa.neut.lymp.mono.eosi.baso.Block.Blockorder/Model2_result.RData")
# sel_criteria(result)
## run the model selection #####
pm1 <- proc.time()
out_pairedend_cbc <-  data.frame(Date=as.Date(character()),
                                   File=character(), 
                                   User=character(), 
                                   stringsAsFactors=FALSE) 

for(criteria in 1){ # i <- 1 criteria <- 1
  model_th <- 1
  full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca +  RINa + 
                               neut + lymp + mono + eosi + baso + 
                               Block + Blockorder)
  repeat{
    out_model <- fit_model(full_model, model_th, criteria)
    #out_model <- s[[i]]
    assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
    ms_val <- get(paste("ms_criteria", model_th, sep = "_" ))
    cov_del <- ms_val[1,criteria] # cov_del <- 14; i <- 1# model_th <- 2 # criteria <- 1
    
    cov_set <- list_model(full_model)$test.mat # dim(cov_set)
    res <- data.frame(criteria = colnames(ms_val)[criteria], 
                      model = model_th, 
                      cov_del = rownames(cov_set)[ cov_del])
    out_pairedend_cbc <- rbind(out_pairedend_cbc, res)
    if (cov_del ==1) break
    block_ind <- grep("Block2", colnames(full_model))
    blockorder_ind <-grep("Blockorder", colnames(full_model))
    indicator <- FALSE
    if(length(block_ind)!=0){
      if(cov_del+1 == cov_set["Block",2])
        {
        full_model <- full_model[, -c(block_ind, block_ind+1, block_ind+2)]
        indicator <- TRUE
      }
    }
    
    if(length(blockorder_ind)!=0){
      if(cov_del+1 == cov_set["Blockorder", 2])
        {
        full_model <- full_model[, -blockorder_ind] 
        indicator <- TRUE
      }
    }
    
    if (indicator == FALSE) full_model <- full_model[, -(cov_del +1)] 
    model_th <- model_th +1
  }
}

# write.csv(out_pairedend_cbc, "out_pairedend_cbc.csv", row.names = F)
proc.time() -pm1

# 
# test.mat

# 
# model_th <- 1000 ## Using QuasiSeq1.0.4#####
# model_th <- 7
# full_model <- model.matrix(~Line +Concb + RINa + 
#                              neut + lymp + mono + baso + 
#                              Block)
# full_model <- model.matrix(~Line + Concb +RINa  +
#                              neut + lymp +mono + baso + Block)
# 
#criteria <- 1 model <- 7
# # out_model <- fit_model(full_model, model_th, 1)
# ## rerun above model 777 with new QuasiSeq####
# ## has some error on glm.fit ###
# 
# list_out <- list_model(full_model)
# 
# design.list <- list_out[[1]]
# 
# #design.list <- list_out$design.list
# #test.mat <- list_out$test.mat
# # fit <- QL.fit(counts, design.list, test.mat, # dim(counts)
# #               log.offset = log.offset, print.progress=FALSE,
# #               Model = "NegBin", method = "optim")
# # result<- QL.results(fit, Plot = FALSE)
# # 
# log.offset <- log(apply(counts, 2, quantile, .75))
# fit2 <- QL.fit(counts, design.list, test.mat,  # dim(counts)
#               log.offset = log.offset, print.progress=TRUE,
#               Model = "NegBin", method = "optim")
# result2<- QL.results(fit2, Plot = FALSE)
# (result[[3]])
# str(result)
# str(fit2)
# str(fit)
# 
# sum((result2$Q.values[[3]]<=.15)&(abs(fit2$coef[,2]/log(2))>=1))
# sum(result2$Q.values[[3]][, "Line"]<=.05)
# hist(result2$P.values[[3]][,"Line"], nclass  =100)
# str(result2$P.values[[3]])
# 
# sum((result$Q.values[[3]][,"Line"]<=.15)&(abs(fit$coef[,2]/log(2))>=1))
# sum(result$Q.values[[3]][,"Line"]<=.05)
# 
# 
# 
# hist(result$P.values[[3]][,"Line"], nclass  =100)
# ##calculate sel_criteria for output models#####
# 
# # dir <- "U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/pvalue05/"
# #dir <- "U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/ks/"
# dir <- "U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue053/"
# dir <- "U:/R/RA/Data/RFI-newdata/resultsimulation/pvalue05/"
# f <- list()
# r <- list()
# s <- list()
# q.line <- vector()
# dir.list <- list.files(dir)
# for(i in 1:10){ # i <- 7
#    diri <- paste0(dir, dir.list[grep(paste0("Model", i,".Line" ), dir.list)], "/")
#   load(file = paste0(diri, paste0("Model", i, "_fit.RData")))
#   load(file = paste0(diri, paste0("Model", i, "_result.RData")))
# #   str(fit)
#   s[[i]] <- sel_criteria(result)
# #   f[[i]] <- fit
# #   r[[i]] <- result
# #   q.line[i] <- sum(r[[i]]$Q.values[[3]][,"Line"]<=0.05)
# #   
# }
# colnames(fit$LRT)
# s[[6]]
# str(fit)
# str(fit2)
# sel_criteria(results)
# str(result$P.values[[3]][,1])
# str(result2$P.values[[3]])
# plot(result$P.values[[3]][,1],result2$P.values[[3]])
# full_model <- model.matrix(~Line)
# model_th <- 10
# 
# listout <- list_model(full_model)
# design.list <- listout$design.list
# test.mat <- listout$test.mat
# test.mat
# 
# fitl <- QL.fit(counts, log.offset = log.offset, Model = "NegBin", 
#                Method = "optim", design.list = design.list, test.mat = test.mat)
# resultl <- QL.results(fitl)
# dev.off()
# hist(resultl$P.values[[3]][,"Line"], nclass = 100)
# sum(resultl$Q.values[[3]][,"Line"] <=.05)
# 
# result <- resultl
# fit <- fitl
# save(result, file = "Model10_result.RData")
# save(fit, file = "Model10_fit.RData")
# 
# load("Model7_result.RData")
# sum(result$Q.values[[3]][,"Line"] <=.05)
# hist(result$P.values[[3]][,"Line"], nclass = 100)
# 
# load("Model1_result.RData")
# 
# sum(result$Q.values[[3]][,"Line"] <=.05)
# hist(result$P.values[[3]][,"Line"], nclass = 100)

