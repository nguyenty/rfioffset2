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
scount <- read.table("paired end uniquely mapped reads count table.txt", 
                     header = T)
scount <- scount[-c(which(scount[,1] %in%"ENSSSCG00000007978"),
                    which(scount[,1] %in%"ENSSSCG00000014725")),]
counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                             rowMeans(scount[,-1])>8 ,-1])

counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                              rowMeans(scount[,-1])>20 ,-1])

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

model_th <- 1
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                             neut + lymp + mono + eosi + baso + 
                             Block + Blockorder)
x1 <- (full_model[,-c(1:2)])
x1s <- scale(x1)
u <- svd(x1s)$u
str(svd(x1s))
cumsum(svd(x1s)$d^2)/sum(svd(x1s)$d^2)

full_model2 <- cbind(full_model[,c(1,2)], u)
design.list <- vector("list", 2)
design.list[[1]] <- full_model2
design.list[[2]] <- full_model2[,-2]


log.offset <- log(apply(counts, 2, quantile, .75))
fittrancov <- QL.fit(counts, design.list = design.list, 
               log.offset = log.offset, 
               print.progress=TRUE,
               Model = "NegBin", method = "optim")
load("Model1_fit.RData")

hist(fittrancov$coef[, 23], nclass = 100, prob = T)
boxplot(fittrancov$coef[, c(2, 3, 23)])
#hist(fit$coef[, 3], nclass = 100, prob = T)
resulttrancov <- QL.results(fittrancov, Plot = FALSE)
hist(fittrancov$coef[,6], nclass = 100, prob = T)
hist(resulttrancov$P.values[[3]], nclass = 100)
pp <- resulttrancov$P.values[[3]]
save(fittrancov, file = "fittrancov.RData")
save(resulttrancov, file = "resulttrancov.RData")

load(file = "fittrancov.RData")
load("resulttrancov.RData")

f <- function(m) class(try(solve(m),silent=T))=="matrix"
se.b1 <- function(b,w,dl2, counts){
  # b is fitted coefficient, w is the negative binomial dispersion
  # dl2 is the design  matrix full model 
  # counts is the counts data
  log.offset <- log(apply(counts, 2, quantile, 0.75))
  X <- dl2
  
  eta <- c(log.offset + X%*%b)
  
  mu <- exp(eta)
  
  fish <- t(X)%*%(diag(nrow(X))*mu/(mu*w+1))%*%X
  
  if(f(fish)){
    fish.inv <- solve(fish)
    return(sqrt(diag(fish.inv)))
  }else{
    return(rep(0, ncol(X)))
  }
}
prefit <- fittrancov
model_th <- 1




ses <- prefit$coefficients

for(j in 1:nrow(prefit$coefficients)){ 
  se.new <- se.b1(prefit$coefficients[j,],prefit$NB.disp[j],full_model2)
  ses[j,] <- se.new
}
dim(ses)

del.gene <- which(apply(ses, 1, sum)==0)
del.gene #  6269  7299  9257 10405 11574




source("quasiseq shrinkage functions2.R")



shrink.phi<-function(phi.hat,den.df){
  phi.hat[phi.hat<=0]<-min(phi.hat[phi.hat>0])
  z<-log(phi.hat); z[z==Inf]<-max(z[z!=Inf]); z[z==-Inf]<-min(z[z!=-Inf]);mnz<-mean(z)
  
  ## solve for d0 and phi0
  d0arg<-var(z)-trigamma((den.df)/2)
  if(d0arg>0){
    dif<-function(x,y) abs(trigamma(x)-y)
    inverse.trigamma<-function(y) optimize(dif,interval=c(0,10000),y=y)$minimum
    d0<-2*inverse.trigamma(d0arg)
    phi0<-exp(mnz-digamma((den.df)/2)+digamma(d0/2)- log(d0/(den.df)))
    
    ## compute shrunken phi's
    phi.shrink<-((den.df)*phi.hat+d0*phi0)/(den.df+d0)  
  }
  else{phi.shrink<-rep(exp(mnz),length(z)); d0<-Inf; phi0<-exp(mnz) }
  return(list(phi.shrink=phi.shrink,d0=d0,phi0=phi0))  
}
#### Code for implementing Smyth's approach ends here ####
load("Model1_fit.RData")
i <- 21
hist(fittrancov$coef[,i], nclass = 100, prob = T)
hist(fit$coef[,i], nclass = 100, prob = T)
hist(ses[,4], nclass = 100)

mu1 <- 0
sigma1 <- .01^2

lf <- function(par,index){
  l <- NULL
  ys <- fittrancov$coef[,index]
  sh <- ses[,index]^2 #index <- 2, dim(ses)
  out <- shrink.phi(sh, 8)[[1]]
  for (i in 1:dim(fittrancov$coef)[1])
    l[i] <- par[1]*dnorm(ys[i], mu1, sqrt(sigma1+out[i])) + 
    (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + out[i]))
  -sum(log(l))
}

par <- c(.5, .1, 1)
outpar <- optim(par = par, lf, lower = c(0, -Inf, 0), upper = c(1, Inf, Inf), index = 2)
outpar$par
?optim
hist(ses[,3], nclass = 100)
####write a function for all this

post.mu <- function(index){
  py <- emu <- NULL
  ys <- fittrancov$coef[,index]
  sh <- ses[,index]^2 #index <- 2, dim(ses)
  out <- shrink.phi(sh, 8)[[1]]
  par <- c(.5, .1, 1)
  outpar <- optim(par = par, lf, lower = c(0, -Inf, 0), upper = c(1, Inf, Inf), index = 2)
  par <- outpar$par
  
  for (i in 1:dim(fittrancov$coef)[1])
  {py[i] <- par[1]*dnorm(ys[i], mu1, sqrt(sigma1+out[i])) + 
     (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + out[i]))
   
   emu[i] <- (par[1]*dnorm(ys[i], mu1, sqrt(sigma1+out[i]))* 
                (1/out[i]*ys[i]+1/sigma1*mu1)/(1/out[i]+1/sigma1)+ 
                (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + out[i]))* 
                (1/out[i]*ys[i]+1/par[3]*par[2])/(1/out[i]+1/par[3]))/py[i]
  }
  emu
}

newbeta <- fittrancov$coef
for(i in 1:dim(newbeta)[2]){
  newbeta[,i] <- post.mu(i)
} 
# hist(newbeta[,9], nclass = 100, prob = T)
# ## prob = 0 is too big####
# 
# 
# 
# #### example beta_3####
# dim(fittrancov$coef)
# y <- fit$coef[, 3]
# sh <- ses[,3]^2
# length(y)
# length(sh)
# mu1 <- 0; sigma1 <- .01^2
# 
# lf <- function(par){
#   l <- NULL
#   for (i in 1:length(y))
#     l[i] <- par[1]*dnorm(y[i], mu1, sqrt(sigma1+sh[i])) + 
#     (1-par[1])*dnorm(y[i], par[2], sqrt(par[3] + sh[i]))
#   -sum(log(l))
# }
# mean(y)
# hist(y, nclass = 100, prob = T)
# 
# pm <- proc.time()
# par <- c(.1, .03, 1)
# outpar <- optim(par = par, lf, lower = c(0, -Inf, 0), upper = c(1, Inf, Inf))
# outpar$par
# proc.time() -pm
# 
# ######################################
# mean(y)
# g <- length(y)
# z <- rbinom(g, 1, .98)
# 
# mu <- z* rnorm(g, 0, 0.000001) + (1-z)*rnorm(g, .2,.000001)
# 
# delta <- sqrt(sh)
# y1 <- rnorm(g, mu, delta)
# hist(y1, nclass = 100)
# hist(y, nclass = 100, prob = T)
# mean(y[-which(y>2)])
# var(y[-which(y>2)])
# curve(dnorm(x, mean(y)/1.3, sd(y)/1.6), add = T, col = "red")
# lines()
# mean(y)
# var(y)
# shapiro.test(y[1:5000])
# boxplot(y)
# plot(log(apply(counts, 1, mean)), sqrt(sh))
# hist(sh, nclass = 100 )
# dim(counts)
# length(sh)
# lf1 <- function(par){
#   l <- NULL
#   for (i in 1:length(y1))
#     l[i] <- par[1]*dnorm(y1[i], 0, sqrt(0.01^2 +sh[i])) + 
#     (1-par[1])*dnorm(y1[i], par[2], sqrt(par[3] + sh[i]))
#   -sum(log(l))
# }
# 
# hist(y, nclass = 100)
# hist(y1, nclass = 100)
# 
# pm <- proc.time()
# par <- c(.9, 2, 1)
# optim(par = c(.1, 2, 1), lf1, lower = c(0, -Inf, 0), upper = c(1, Inf, Inf))
# proc.time()-pm
# 
# 
# load("new.beta.RData")
# str(new.beta)
# 
# bvec <- matrix(0, nrow = length(new.beta[[1]]), ncol = length(new.beta))
# 
# for(i in 1:length(new.beta)){
#   bvec[,i] <- new.beta[[i]]
# }
# 
# model_th <- 11111 # test for Emperial Bayes method#####
# full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
#                              lneut + llymp + lmono + leosi + lbaso + 
#                              Block + Blockorder)
# new.offset <- t(full_model[, -2]%*% t(bvec[,-2]))
# dim(new.offset)
# 
# plot(newmu, fit$coef[,2])
# 
# 
# y <- fit$coef[,2]
# phi.hat <- ses[,2]^2
# den.df <- 8
# out <- shrink.phi(phi.hat, den.df)
# 
# ################
# #hist(y,nclass = 100)
# library("plyr")
# mu1 <- 0
# sigma1 <- .001^2
# 
# lf <- function(par){
#   l <- NULL
#   for (i in 1:length(y))
#     l[i] <- par[1]*dnorm(y[i], mu1, sqrt(sigma1+out[[1]][i])) + 
#     (1-par[1])*dnorm(y[i], par[2], sqrt(par[3] + out[[1]][i]))
#   -sum(log(l))
# }
# 
# # mean(y)
# # hist(y, nclass = 100, prob = T)
# 
# pm <- proc.time()
# par <- c(.5, .5, 1)
# outpar <- optim(par = par, lf, lower = c(0, -Inf, 0), upper = c(1, Inf, Inf))
# outpar$par
# proc.time() -pm
# 
# 
# ######simulation#####
# 
# phi.hat <- ses[,2]^2
# den.df <- 8
# out <- shrink.phi(phi.hat, den.df)
# 
# z <-rbinom(length(y), 1, .9)
# mu11 <- 0
# sigma11 <- .01^2
# sigma12 <- 1
# mu12 <- 2
# ys  <- z*rnorm(length(y), mu11, sqrt(sigma11+out[[1]])) +
#   (1-z)*rnorm(length(y), mu12, sqrt(sigma12+out[[1]]))
# 
# sigma1 <- .01^2
# 
# 
# lfs <- function(par){
#   l <- NULL
#   for (i in 1:length(y))
#     l[i] <- par[1]*dnorm(ys[i], mu11, sqrt(sigma11+out[[1]][i])) + 
#     (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + out[[1]][i]))
#   -sum(log(l))
# }
# 
# pm <- proc.time()
# par <- c(.9, .1, .5)
# outpar <- optim(par = par, lfs, lower = c(0, -Inf, 0), upper = c(1, Inf, Inf))
# outpar$par
# proc.time() -pm
# 
# lf <- function(par,index){
#   l <- NULL
#   ys <- fit$coef[,index]
#   sh <- ses[,index]^2 #index <- 2, dim(ses)
#   out <- shrink.phi(sh, 8)[[1]]
#   for (i in 1:dim(fit$coef))
#     l[i] <- par[1]*dnorm(ys[i], mu1, sqrt(sigma1+out[i])) + 
#     (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + out[i]))
#   -sum(log(l))
# }
# 
# par <- c(.5, .5, 1)
# outpar <- optim(par = par, lf, lower = c(0, -Inf, 0), upper = c(1, Inf, Inf), index = 2)
# outpar$par
# 
# ####write a function for all this
# 
# post.mu <- function(par,index){
#   py <- emu <- NULL
#   ys <- fit$coef[,index]
#   sh <- ses[,index]^2 #index <- 2, dim(ses)
#   sh2 <- shrink.phi(sh, 8)[[1]]
#   for (i in 1:dim(fit$coef))
#     {py[i] <- par[1]*dnorm(ys[i], mu1, sqrt(sigma1+out[[1]][i])) + 
#     (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + out[[1]][i]))
#   
#   emu[i] <- (par[1]*dnorm(ys[i], mu1, sqrt(sigma1+out[[1]][i]))* (1/sh*ys[i]+1/sigma1*mu1)/(1/sh+1/sigma1)+ 
#     (1-par[1])*dnorm(ys[i], par[2], sqrt(par[3] + out[[1]][i]))* (1/sh*ys[i]+1/par[3]*par[2])/(1/sh+1/par[3]))/y[i]
# }
# emu
# }
# 
# newmu <- post.mu(outpar$par, 2)
# plot(newmu, fit$coef[,2])
# ### estimate of pi is 1 for all the case.#####




### shrink ses######
newses <- ses
for (i in 1:dim(ses)[2]){
  phi.hat <- ses[,i]^2
  den.df <- 8
  out <- shrink.phi(phi.hat, den.df)
  newses[,i] <- out[[1]]
}

## Using bayesian with horseshoe prior, unshrinked ses #######
library(rjags)
#########modelh##############
modelh <- "
model{
# likelihood 
for (i in 1:length(y)){
y[i] ~ dnorm(mu[i], 1/sigma[i])
mu[i] ~ dnorm(0, 1/delta[i])
delta[i] ~dt(0, 1/s^2, 1) T(0,)

}
# prior level 1
s ~ dt(0, 5, 1) T(0,)
}
"
dim(fittrancov$coef)
newcoef2 <- fittrancov$coef[,-c(1,2)]
pm <- proc.time() 
for (i in 3:dim(newses)[2])
{
  data <- list(y =  fittrancov$coef[,i], 
             sigma = newses[,i]^2)
  mm <- jags.model(textConnection(modelh), data, n.adapt=100,n.chains = 1) # mix point mass
  resm <- coda.samples(mm, c("mu"), 500) # mix point mass
  pos <- apply(resm[[1]], 2, mean)
  newcoef2[,i-2] <- pos
}
proc.time() -pm
save(newcoef2, file = "newcoef2.RData")
load("newcoef2.RData")
str(newcoef2)
bayesbeta <- newcoef2%*% t(full_model2[, -c(1,2)])
newcoef3 <- newcoef2
head(newcoef2)
head(newcoef3)
bayesbeta3 <- newcoef3%*% t(full_model2[, -c(1,2)])

design.list.beta <- vector("list", 2)
design.list.beta[[1]] <- model.matrix(~Line) # line
design.list.beta[[2]] <-rep(1, 31)
fitbayesbeta <- myQL.fit(counts, design.list = design.list.beta, 
                         betavec = bayesbeta3,   # dim(counts)
                    log.offset = log.offset, 
                    print.progress=TRUE,
                    Model = "NegBin", method = "optim")

save(fitbayesbeta, file = "fitbayesbeta.RData")
resultbayesbeta <- QL.results(fitbayesbeta, Plot = FALSE)
save(resultbayesbeta, file = "resultbayesbeta.RData")
load(file = "resultbayesbeta.RData")
hist(resultbayesbeta$P.values[[3]], nclass = 100)
sum(resultbayesbeta$Q.values[[3]]<=.10)

str(resm)
dim(newcoef2)
plot(fittrancov$coef[,4],newcoef2[,2])
# hist(pos, nclass = 1000)

mh <- jags.model(textConnection(model1), data, n.chains = 1) # horseshoe

resh <- coda.samples(mh, c("tau", "beta"), 5000) # horseshoe
proc.time() -m0

######fit myQLfit with a matrix of additional offset, do not shrink coefficients####
source("quasiseq shrinkage functions2.R")

betavec <- fittrancov$coef[,-c(1,2)]%*%t(full_model2[,-c(1,2)])
dim(betavec)

design.list.beta <- vector("list", 2)
design.list.beta[[1]] <- model.matrix(~Line) # line
design.list.beta[[2]] <-rep(1, 31)

fitbeta <- myQL.fit(counts, design.list = design.list.beta, betavec = betavec,   # dim(counts)
              log.offset = log.offset, 
              print.progress=TRUE,
              Model = "NegBin", method = "optim")
str(fitbeta)
save(fitbeta, file = "fitbeta.RData")
resultbeta <- QL.results(fitbeta, Plot = FALSE)
save(resultbeta, file = "resultbeta.RData")

hist(resultbeta$P.values[[3]], nclass = 100)
sum(resultbeta$Q.values[[3]]<=.01)



########run qlfit with best model from model selection #####

full_model <- model.matrix(~Line + RINa + Concb +neut + lymp + mono + Block)
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
log.offset <- log(apply(counts, 2, quantile, .75))
fitbest <- QL.fit(counts, design.list, test.mat, # dim(counts)
              log.offset = log.offset, 
              print.progress=TRUE,
              Model = "NegBin", method = "optim")

resultbest <- QL.results(fitbest, Plot = FALSE)
save(fitbest, file = "fitbest.RData")
save(resultbest, file = "resultbest.RData")
load(file = "fitbest.RData")
load(file = "resultbest.RData")
dim(resultbest$P.values[[3]])

hist(resultbest$P.values[[3]][,"Line"], nclass = 100)
sum(resultbest$Q.values[[3]][,"Line"] <= .05)
str(result)
hist(result$P.values[[3]][, "Line"], nclass =100)
plot(-log(resultbest$P.values[[3]][, "Line"]), -log(resultnewbeta$P.values[[3]]))


