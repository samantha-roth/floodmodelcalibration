
#Calibration

rm(list=ls())
library(mvtnorm); library(invgamma); library(spam)
#source functions from Kennedy et al 2021
source("C:/FloodingModelCalibrationProject/sml-athena-main/GPfunctionsOptim.R")
source("C:/FloodingModelCalibrationProject/sml-athena-main/hetGPfunctions.R")

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/pars.homgp.cheap.RData")


## define test metrics
MSE <- function(true, pred){
  mean((true-pred)^2)
}

Score <- function(y, m, v){
  # y = "true" value
  # m = emulator mean
  # v = emulator variance - including nugget term
  -((y-m)^2)/v - log(v)
}


#load predictions
parsTrue<- read.csv("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/parVals/RunTrue_1.csv")
parsTrue<- as.data.frame(parsTrue)

#load parameters
load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/allParVals.RData")
parVals10m<- as.data.frame(parVals)
load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/allParVals.RData")
parVals50m<- as.data.frame(parVals)


ch.true<- parsTrue$n_ch
fp.true<- parsTrue$n_fp
rwe.true<- parsTrue$rwe
ree.true<- parsTrue$ree

#load predictions
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10m.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10mfrom50m.RData")

nPars=4

##training data

## load data
#y data
y.test1 <- 0 #true Euclidean distance between output and reality
y.gp1<- metrics.10mfrom50m$EuclideanDists
y.c1<- metrics.10mfrom50m$EuclideanDists
y.e1<- metrics.10m$EuclideanDists

y.test <- y.test1
y.gp <- y.gp1
y.c <- y.c1
y.e <- y.e1

z=y.test

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2)
n.e <- length(y.e2)
n.ml <- length(y.e)

#x data
x.true<- as.matrix(parsTrue[,-(1:2)])

x.gp<- as.matrix(parVals50m[,-1])
x.c <- as.matrix(parVals50m[,-1])
x.e <- as.matrix(parVals10m[,-1])

## standardise inputs (x-mu)/sig
x.e.std <- scale(x.e)
x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
x.std <- scale(x.gp)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"
x.c.std.rev <- x.c.std[n.c:1,]
x.e.std.rev <- x.e.std[n.ml:1,]

##get example data for illustrating the fit
#library(lhs)
#nTest<- 20
#set.seed(14)
#lh.Test<- lhs::maximinLHS(nTest,nPars)
##A 100 by 2 Latin Hypercube Sample matrix with values uniformly distributed on [0,1]
##transform into samples from the right uniform distributions
#samp.Test<- matrix(0, nTest, nPars)
#samp.Test[,1]<- qunif(lh.Test[,1], min = 0.02, max = 0.1)
#samp.Test[,2]<- qunif(lh.Test[,2], min = 0.02, max = 0.4)
##samp.Test[,3]<- qunif(lh.Test[,3], min = -1, max= 1)
##samp.Test[,3]<- qbeta(lh.Test[,3], 10, 10)
#samp.Test[,3]<- qnsbeta(lh.Test[,3], 10, 10, min = 0, max = 2)
#samp.Test[,4]<- qnorm(lh.Test[,4], mean = 0, sd = 3)

#x.test<- samp.Test
x.test<- x.true

#quantities for homGP emulator

var.beta <- diag(1, ncol(x.test)+1)
mean.beta <- rep(0, ncol(x.test)+1)

H <- cbind(1, x.gp)

covmat.gp <- H %*% var.beta %*% t(H) + 
  cov.x1.x2(x.std, x.std, pars$m_sigma, pars$m_theta, 2) + diag(rep(pars$m_nugget),length(y.gp))

precmat.gp <- chol2inv(chol(covmat.gp))
precYminusHb<- (precmat.gp %*% (y.gp - H %*% mean.beta))


################################################################################
#gibbs sampler for sigma^2

#Code to generate a random sample from the posterior distribution of sigma2
rpost_sigma2<- function(mean.ml,z){
  prior<-c(0.2,0.2)
  sh= nrow(x.test)/2 + prior[1]
  ra= .5*(z-mean.ml)^2 + prior[2]
  invgamma::rinvgamma(1,shape=sh,rate=ra) 
}
################################################################################

#range= c(ch.true-.001,ch.true+.001)
#colnames(metrics.10m)[4:6]<- c("EuclideanDists","Fvals","Cvals")
#colnames(metrics.10mfrom50m)[4:6]<- c("EuclideanDists","Fvals","Cvals")
#allmetrics<- rbind(metrics.10m[-test.inds,],metrics.10mfrom50m[-test.inds,])
#close.inds<- intersect(which(allmetrics$n_ch> range[1]),which(allmetrics$n_ch< range[2]))
#fake.delta<- mean(allmetrics$EuclideanDists[close.inds])

#set hyperparameters
#sigma2REE= 3
#alpha_RWE= 10
#beta_RWE=10

#run the code for the Gibbs sampler#Code the Gibbs sampler.
mh.alg<- function(init, n.sample) {
  x1.t <- init[1] #ch
  x2.t <- init[2] #fp
  x3.t <- init[3] #ch
  x4.t <- init[4] #fp
  x5.t <- init[5] #sigma^2
  
  x1.out <- rep(NA,n.sample+1)
  x2.out <- rep(NA,n.sample+1)
  x3.out <- rep(NA,n.sample+1)
  x4.out <- rep(NA,n.sample+1)
  x5.out <- rep(NA,n.sample+1)
  
  x1.out[1]= x1.t
  x2.out[1]= x2.t
  x3.out[1]= x3.t
  x4.out[1]= x4.t
  x5.out[1]= x5.t
  
  for (i in 1 : n.sample) {
    #print(paste("Iteration",i))
    if(i%%1000==0){print(paste("Iteration",i))}
    #sample from the posterior distributions
    
    x.test<- matrix(c(x1.out[i], x2.out[i],x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i], x2.out[i],x3.out[i], x4.out[i])))
    
    ##scale validation data accordingly
    x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
    
    #calculated quantities after fitting homSML emulator
    H.new <- cbind(1, x.test)
    crosscov <- gp.crosscov(x.v1, x.std, pars$m_sigma, pars$m_theta, H.new, H, var.beta )
    
    mean.gp <- crosscov %*% precYminusHb
    
    
    #MH algorithm for ch
    u<- runif(1)
    #prop<- rnorm(1,mean=.02+ (.1-.02)/2,sd= (.1-.02)/2) #symmetric proposal
    prop<- runif(1,min=.02,max=.1)
    if(prop<.02 | prop>.1){
      x1.out[i+1]<- x1.out[i]
    } else{
      #compute the inverse covariance matrix of the deltas
      
      prop.x.test<- matrix(c(prop, x2.out[i], x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i], x2.out[i], x3.out[i], x4.out[i])))
      
      prop.x.v1 <- scale(prop.x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
      
      #calculated quantities after fitting homSML emulator
      prop.H.new <- cbind(1, prop.x.test)
      prop.crosscov <- gp.crosscov(prop.x.v1, x.std, pars$m_sigma, pars$m_theta, prop.H.new, H, var.beta )
      
      prop.mean.gp <- prop.crosscov %*% precYminusHb
      
      
      logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.gp)^2
      logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.gp)^2
      
      logMetRatio<- logMetRatio1-logMetRatio2
      
      #metropolis ratio in the metropolis hastings algorithm 
      if(log(u)<=logMetRatio) x1.out[i+1]<- prop
      if(log(u)>logMetRatio) x1.out[i+1]<- x1.out[i]
    }
    
    #if() x1.out[i+1]<- x1.out[i]
    #if(prop>=.02 && prop<=.1){
    
    x.test<- matrix(c(x1.out[i+1], x2.out[i], x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i+1], x2.out[i],x3.out[i], x4.out[i])))
    
    ##scale validation data accordingly
    x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
    
    #calculated quantities after fitting homSML emulator
    H.new <- cbind(1, x.test)
    crosscov <- gp.crosscov(x.v1, x.std, pars$m_sigma, pars$m_theta, H.new, H, var.beta )
    
    mean.gp <- crosscov %*% precYminusHb
    
    #MH algorithm for fp
    u<- runif(1)
    #prop<- rnorm(1,mean=.02+((4.-.02)/2),sd= (.4-.02)/2) #symmetric proposal
    prop<- runif(1,min=.02,max=.4)
    if(prop<.02 | prop>.4){
      x2.out[i+1]<- x2.out[i]
    } else{
      prop.x.test<- matrix(c(x1.out[i+1], prop, x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i], x2.out[i], x3.out[i], x4.out[i])))
      
      prop.x.v1 <- scale(prop.x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
      
      #calculated quantities after fitting homSML emulator
      prop.H.new <- cbind(1, prop.x.test)
      prop.crosscov <- gp.crosscov(prop.x.v1, x.std, pars$m_sigma, pars$m_theta, prop.H.new, H, var.beta )
      
      prop.mean.gp <- prop.crosscov %*% precYminusHb
      
      logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.gp)^2
      logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.gp)^2
      
      logMetRatio<- logMetRatio1-logMetRatio2
      
      #metropolis ratio in the metropolis hastings algorithm 
      if(log(u)<=logMetRatio) x2.out[i+1]<- prop
      if(log(u)>logMetRatio) x2.out[i+1]<- x2.out[i]
    }
    #if() x2.out[i+1]<- x2.out[i]
    #if(prop>=.02 && prop<= .4){
    
    x.test<- matrix(c(x1.out[i+1], x2.out[i+1], x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i+1], x2.out[i+1], x3.out[i], x4.out[i])))
    
    ##scale validation data accordingly
    x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
    
    #calculated quantities after fitting homSML emulator
    H.new <- cbind(1, x.test)
    crosscov <- gp.crosscov(x.v1, x.std, pars$m_sigma, pars$m_theta, H.new, H, var.beta )
    
    mean.gp <- crosscov %*% precYminusHb
    
    
    #MH algorithm for RWE
    u<- runif(1)
    #prop<- rnorm(1,mean=1,sd= 1) #symmetric proposal
    prop<- runif(1,min=.95,max=1.05)
    if(prop<.95 | prop>1.05){
      x3.out[i+1]<- x3.out[i]
    } else{
      prop.x.test<- matrix(c(x1.out[i+1],x2.out[i+1], prop, x4.out[i]),nrow=1,ncol= length(c(x1.out[i+1],x2.out[i+1], prop, x4.out[i])))
      
      prop.x.v1 <- scale(prop.x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
      
      #calculated quantities after fitting homSML emulator
      prop.H.new <- cbind(1, prop.x.test)
      prop.crosscov <- gp.crosscov(prop.x.v1, x.std, pars$m_sigma, pars$m_theta, prop.H.new, H, var.beta )
      
      prop.mean.gp <- prop.crosscov %*% precYminusHb
      
      #logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.gp)^2 + (alpha_RWE-1)*log(prop) + (beta_RWE-1)*log(2-prop)
      #logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.gp)^2 + (alpha_RWE-1)*log(x3.out[i]) + (beta_RWE-1)*log(2-x3.out[i])
      
      logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.gp)^2
      logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.gp)^2 
      
      logMetRatio<- logMetRatio1-logMetRatio2
      
      #metropolis ratio in the metropolis hastings algorithm 
      if(log(u)<=logMetRatio) x3.out[i+1]<- prop
      if(log(u)>logMetRatio) x3.out[i+1]<- x3.out[i]
    }
    
    x.test<- matrix(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], x4.out[i]),nrow=1,ncol= length(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], x4.out[i])))
    
    ##scale validation data accordingly
    x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
    
    #calculated quantities after fitting homSML emulator
    H.new <- cbind(1, x.test)
    crosscov <- gp.crosscov(x.v1, x.std, pars$m_sigma, pars$m_theta, H.new, H, var.beta )
    
    mean.gp <- crosscov %*% precYminusHb
    
    #MH algorithm for REE
    u<- runif(1)
    #prop<- rnorm(1,mean=0,sd= 3) #symmetric proposal
    prop<- runif(1,min=-5,max=5)
    #if(prop<0 | prop>2){
    #  x4.out[i+1]<- x4.out[i]
    #} else{
    prop.x.test<- matrix(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], prop),nrow=1,ncol= length(c(x1.out[i+1],x2.out[i+1], x3.out[i+1], prop)))
    
    prop.x.v1 <- scale(prop.x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
    
    #calculated quantities after fitting homSML emulator
    prop.H.new <- cbind(1, prop.x.test)
    prop.crosscov <- gp.crosscov(prop.x.v1, x.std, pars$m_sigma, pars$m_theta, prop.H.new, H, var.beta )
    
    prop.mean.gp <- prop.crosscov %*% precYminusHb
    
    #logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.gp)^2 - (1/(2*sigma2REE))*prop^2
    #logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.gp)^2 - (1/(2*sigma2REE))*(x4.out[i])^2
    
    logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.gp)^2 
    logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.gp)^2 
    
    logMetRatio<- logMetRatio1-logMetRatio2
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio) x4.out[i+1]<- prop
    if(log(u)>logMetRatio) x4.out[i+1]<- x4.out[i]
    #}
    
    x.test<- matrix(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], x4.out[i+1]),nrow=1,ncol= length(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], x4.out[i+1])))
    
    ##scale validation data accordingly
    x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
    
    #calculated quantities after fitting homSML emulator
    H.new <- cbind(1, x.test)
    crosscov <- gp.crosscov(x.v1, x.std, pars$m_sigma, pars$m_theta, H.new, H, var.beta )
    
    mean.gp <- crosscov %*% precYminusHb
    
    
    #sigma^2
    x5.out[i+1] <-rpost_sigma2(mean.ml=mean.gp, z=z)
  }
  out <- cbind(x1.out, x2.out, x3.out, x4.out, x5.out)
  out
}


#set.seed(142)
#init_vals<- c(runif(1,.02,.1),runif(1,.02,.4),invgamma::rinvgamma(1,.2,.2))
#init_vals<- c(.02+ (.1-.02)/2, .02+ (.4-.02)/2, 1, 0, .2/(.2+1)) #use mode of IG dist

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap.RData")
init_vals<-  as.numeric(res[nrow(res),]); rm(res)


pt<-proc.time() # Start Time
#res <- mh.alg(init = init_vals, n.sample = 1e5)
res <- mh.alg(init = init_vals, n.sample = 2e5)
ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3] # End Time to be used in Effective Samples per Second Calculation

#save(res,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap.RData")
#save(ptFinal,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/time_homGP_cheap.RData")

save(res,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap2.RData")
save(ptFinal,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/time_homGP_cheap2.RData")


load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap.RData")
res1<- res
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap2.RData")
res2<- res

res<- rbind(res1,res2)

nPars=5
for(i in 1:nPars){
  plot(density(res[,i]),main=paste("HomGP Cheap variable ",i,sep=""),col="black")
  lines(density(res[1:2e5,i]),col="red")
}

save(res,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap.RData")

#EDIT FROM HERE
################################################################################
################################################################################
####################Performance evaluation######################################
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP.RData")

hist(res[,1],main="Ch Density",xlab="Ch Values")
# Changing x axis
xtick<-seq(0.02, 0.1, by=0.01)
axis(side=1, at=xtick, labels = FALSE)
hist(res[,2],main="Fp Density",xlab="Fp Values")
hist(res[,3],main="RWE Density",xlab="RWE Values")
hist(res[,4],main="REE Density",xlab="REE Values")
hist(res[,5],main="Sigma^2 density",xlab="Sigma^2 Values")

plot(1:length(res[,1]),res[,1],xlab="Iteration",ylab="Channel roughness")
plot(1:length(res[,2]),res[,2],xlab="Iteration",ylab="Floodplain roughness")
plot(1:length(res[,3]),res[,3],xlab="Iteration",ylab="River width error")
plot(1:length(res[,3]),res[,4],xlab="Iteration",ylab="Riverbed elevation error")
plot(1:length(res[,5]),res[,5],xlab="Iteration",ylab="Sigma^2")


mean(res[,1]); median(res[,1])
mean(res[,2]); median(res[,2])
mean(res[,3]); median(res[,3])
mean(res[,4]); median(res[,4])
mean(res[,5]); median(res[,5])


################################################################################
#get predictions from the emulator at the fit parameter values

fitPars<- apply(res[,1:4],2,median)

fit.x.test<- matrix(fitPars,nrow=1,ncol= length(fitPars))

##scale validation data accordingly
fit.x.v1 <- scale(fit.x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))

#calculated quantities after fitting homSML emulator
fit.H.new <- cbind(1, fit.x.test)
fit.crosscov <- gp.crosscov(fit.x.v1, x.std, pars$m_sigma, pars$m_theta, fit.H.new, H, var.beta )

fit.mean.gp <- fit.crosscov %*% precYminusHb

mc.mean.gp<- rep(NA,nrow(res))
for(i in 1:nrow(res)){
  
  mc.x.test<- matrix(res[i,1:4],nrow=1,ncol= length(res[i,1:4]))
  
  ##scale validation data accordingly
  mc.x.v1 <- scale(mc.x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
  
  #calculated quantities after fitting homSML emulator
  mc.H.new <- cbind(1, mc.x.test)
  mc.crosscov <- gp.crosscov(mc.x.v1, x.std, pars$m_sigma, pars$m_theta, mc.H.new, H, var.beta )
  
  mc.mean.gp[i] <- mc.crosscov %*% precYminusHb
}

hist(mc.mean.gp,breaks=20,main="Histogram of predicted Euclidean distance at each MC sample",xlab="Predicted Euclidean Distance")
xtick<-seq(50, 1400, by=50)
axis(side=1, at=xtick, labels = FALSE)

mean(mc.mean.gp); median(mc.mean.gp)

################################################################################
#distribution of predicted euclidean distances from priors of parameters

load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/parVals/allParVals.RData")
parVals50m<- as.data.frame(parVals)
priorVals<- parVals50m[,-1]

prior.mean.ml<- rep(NA,nrow(priorVals))
for(i in 1:nrow(priorVals)){
  
  prior.x.test<- matrix(priorVals[i,],nrow=1,ncol= length(priorVals[i,]))
  
  ##scale validation data accordingly
  prior.x.v1 <- scale(prior.x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
  
  #calculated quantities after fitting homSML emulator
  prior.H.new <- cbind(1, prior.x.test)
  prior.crosscov <- gp.crosscov(prior.x.v1, x.std, pars$m_sigma, pars$m_theta, prior.H.new, H, var.beta )
  
  prior.mean.gp[i] <- prior.crosscov %*% precYminusHb
}

hist(prior.mean.gp,breaks=20,main="Histogram of predicted Euclidean distance at each prior sample",xlab="Predicted Euclidean Distance")
xtick<-seq(50, 1400, by=50)
axis(side=1, at=xtick, labels = FALSE)

mean(prior.mean.gp); median(prior.mean.gp)
