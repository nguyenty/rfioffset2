# biocLite("edgeR")
library(Matrix)
library(edgeR)
library(reshape)
library(plyr)
library(fields)
library(reshape)
library(fdrtool)
library(QuasiSeq)
# source("QL.fit2.R")
# source("QL.results.R")
# source("NBDev4.R")
# source("PoisDev.R")
# source("negbin.br.R")
# source("fbrglm.R")
# source("glm.fit0.R")
library(AUC)
scount <- read.table("paired end uniquely mapped reads count table.txt", 
                     header = T)
scount <- scount[-c(which(scount[,1] %in%"ENSSSCG00000007978"),
                    which(scount[,1] %in%"ENSSSCG00000014725")),]
# counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
#                              rowMeans(scount[,-1])>8 ,-1])

counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                              rowMeans(scount[,-1])>20 ,-1])

log.offset <- log(apply(counts, 2, quantile, 0.75))
covset <- read.table("covset.txt")
Blockorder <- as.factor(covset$Blockorder)
Block <- as.factor(covset$Block)
Line <- as.factor(covset$Line)
Diet <- as.factor(covset$Diet)
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
resultdir <- getwd()
