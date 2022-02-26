#Calibration- HomML with no delta estimate

rm(list=ls())
library(invgamma)
#source functions from Kennedy et al 2021
source("C:/FloodingModelCalibrationProject/sml-athena-main/GPfunctionsOptim.R")
source("C:/FloodingModelCalibrationProject/sml-athena-main/hetGPfunctions.R")

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/pars.homML.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/pars.homML.RData")


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

logit<- function(x){
  log(x/(1-x))
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
y.gp1<- metrics.10m$EuclideanDists
y.c1<- metrics.10mfrom50m$EuclideanDists
y.e1<- metrics.10m$EuclideanDists

y.test <- y.test1
y.gp <- y.gp1
y.c <- y.c1
y.e <- y.e1

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2); n.e <- length(y.e2)
n.ml <- length(y.e)

#x data
x.true<- as.matrix(parsTrue[,-(1:2)])

x.gp<- as.matrix(parVals10m[,-1])
x.c <- as.matrix(parVals50m[,-1])
x.e <- as.matrix(parVals10m[,-1])

## standardise inputs (x-mu)/sig
x.e.std <- scale(x.e)
x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
x.std <- scale(x.gp)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"
x.c.std.rev <- x.c.std[n.c:1,]
x.e.std.rev <- x.e.std[n.ml:1,]

#get example data for illustrating the fit
#library(lhs); library(extraDistr)
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

#quantities for homSML emulator

m.h <- cbind(cbind(1, x.c[n.c:1,]), rbind(matrix(0, ncol = (ncol(x.e)+1), nrow = n.c - n.e) ,cbind(1, x.e[n.e:1,])))
var.beta <- diag(1, ncol(m.h))

#calculated quantities after fitting homSML emulator
designmat <- matrix(0, ncol=ncol(m.h), nrow=length(c(y.c, y.e)))
designmat[1:length(y.c),1:(ncol(x.e)+1)] <- cbind(1, x.c[n.c:1,])
designmat[-(1:length(y.c)), (1:(ncol(x.e)+1))] <- pars$rho*cbind(1, x.e[n.e:1,])
designmat[-(1:length(y.c)), -(1:(ncol(x.e)+1))] <- cbind(1, x.e[n.e:1,])
dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
precmat.ml <- chol2inv(chol(dat.covmat))
rm(dat.covmat)
H.c <- cbind(1, x.c[n.c:1,])
H.e <- cbind(1, x.c[n.e:1,])

var.bc <- var.beta[1:(ncol(x.e)+1),1:(ncol(x.e)+1)]
var.be <- var.beta[(ncol(x.e)+2):(ncol(m.h)), (ncol(x.e)+2):(ncol(m.h))]

z=y.test

precmat.mlY<- (precmat.ml %*%(c(y.c2, y.e2)))

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
    x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    H <- cbind(1, x.test)
    crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    mean.ml <- crosscov.ml %*% precmat.mlY
    
    #MH algorithm for ch
    u<- runif(1)
    #prop<- rnorm(1,mean= logit((x1.out[i]-.02)/(.1-.02)), sd= 1.7) #symmetric proposal
    prop<- runif(1,min=.02,max=.1)
    #if(prop<.02 | prop>.1){
    #  x1.out[i+1]<- x1.out[i]
    #} else{
    #compute the inverse covariance matrix of the deltas
    
    prop.x.test<- matrix(c(prop, x2.out[i], x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i], x2.out[i], x3.out[i], x4.out[i])))
    
    ##scale validation data accordingly
    prop.x.v2 <- scale(prop.x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    prop.H <- cbind(1, prop.x.test)
    prop.crosscov.ml <- ml.crosscov(prop.x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, prop.H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    prop.mean.ml <- prop.crosscov.ml %*% precmat.mlY
    
    logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.ml)^2
    logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.ml)^2
    
    logMetRatio<- logMetRatio1-logMetRatio2
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio) x1.out[i+1]<- prop
    if(log(u)>logMetRatio) x1.out[i+1]<- x1.out[i]
    #}
    
    #if() x1.out[i+1]<- x1.out[i]
    #if(prop>=.02 && prop<=.1){
    
    x.test<- matrix(c(x1.out[i+1], x2.out[i], x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i+1], x2.out[i],x3.out[i], x4.out[i])))
    
    ##scale validation data accordingly
    x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    H <- cbind(1, x.test)
    crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    mean.ml <- crosscov.ml %*% precmat.mlY
    
    #MH algorithm for fp
    u<- runif(1)
    #prop<- rnorm(1,mean=.02+((.4-.02)/2),sd= (.4-.02)/2) #symmetric proposal
    prop<- runif(1,min=.02,max=.4)
    #if(prop<.02 | prop>.4){
    #  x2.out[i+1]<- x2.out[i]
    #} else{
    prop.x.test<- matrix(c(x1.out[i+1], prop, x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i], x2.out[i], x3.out[i], x4.out[i])))
    
    ##scale validation data accordingly
    prop.x.v2 <- scale(prop.x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    prop.H <- cbind(1, prop.x.test)
    prop.crosscov.ml <- ml.crosscov(prop.x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, prop.H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    prop.mean.ml <- prop.crosscov.ml %*% precmat.mlY
    
    logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.ml)^2
    logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.ml)^2
    
    logMetRatio<- logMetRatio1-logMetRatio2
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio) x2.out[i+1]<- prop
    if(log(u)>logMetRatio) x2.out[i+1]<- x2.out[i]
    #}
    #if() x2.out[i+1]<- x2.out[i]
    #if(prop>=.02 && prop<= .4){
    
    x.test<- matrix(c(x1.out[i+1], x2.out[i+1], x3.out[i], x4.out[i]),nrow=1,ncol= length(c(x1.out[i+1], x2.out[i+1], x3.out[i], x4.out[i])))
    
    ##scale validation data accordingly
    x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    H <- cbind(1, x.test)
    crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    mean.ml <- crosscov.ml %*% precmat.mlY
    
    
    #MH algorithm for RWE
    u<- runif(1)
    #prop<- rnorm(1,mean=1,sd= 1) #symmetric proposal
    prop<- runif(1,min=.95,max=1.05)
    #if(prop<.95 | prop>1.05){
    #  x3.out[i+1]<- x3.out[i]
    #} else{
    prop.x.test<- matrix(c(x1.out[i+1], x2.out[i+1], prop, x4.out[i]),nrow=1,ncol= length(c(x1.out[i+1],x2.out[i+1], prop, x4.out[i])))
    
    ##scale validation data accordingly
    prop.x.v2 <- scale(prop.x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    prop.H <- cbind(1, prop.x.test)
    prop.crosscov.ml <- ml.crosscov(prop.x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, prop.H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    prop.mean.ml <- prop.crosscov.ml %*% precmat.mlY
    
    logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.ml)^2
    logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.ml)^2
    
    logMetRatio<- logMetRatio1-logMetRatio2
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio) x3.out[i+1]<- prop
    if(log(u)>logMetRatio) x3.out[i+1]<- x3.out[i]
    #}
    
    x.test<- matrix(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], x4.out[i]),nrow=1,ncol= length(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], x4.out[i])))
    
    ##scale validation data accordingly
    x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    H <- cbind(1, x.test)
    crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    mean.ml <- crosscov.ml %*% precmat.mlY
    
    #MH algorithm for REE
    u<- runif(1)
    #prop<- rnorm(1,mean=0,sd= 3) #symmetric proposal
    prop<- runif(1,min=-5,max=5)
    #if(prop<-5 | prop>5){
    #  x4.out[i+1]<- x4.out[i]
    #} else{
    prop.x.test<- matrix(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], prop),nrow=1,ncol= length(c(x1.out[i+1],x2.out[i+1], x3.out[i+1], prop)))
    
    ##scale validation data accordingly
    prop.x.v2 <- scale(prop.x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    prop.H <- cbind(1, prop.x.test)
    prop.crosscov.ml <- ml.crosscov(prop.x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, prop.H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    prop.mean.ml <- prop.crosscov.ml %*% precmat.mlY
    
    logMetRatio1<- -(1/(2*x5.out[i]))*(z-prop.mean.ml)^2
    logMetRatio2<- -(1/(2*x5.out[i]))*(z-mean.ml)^2
    
    logMetRatio<- logMetRatio1-logMetRatio2
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio) x4.out[i+1]<- prop
    if(log(u)>logMetRatio) x4.out[i+1]<- x4.out[i]
    #}
    
    x.test<- matrix(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], x4.out[i+1]),nrow=1,ncol= length(c(x1.out[i+1], x2.out[i+1], x3.out[i+1], x4.out[i+1])))
    
    ##scale validation data accordingly
    x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
    
    #calculated quantities after fitting homSML emulator
    H <- cbind(1, x.test)
    crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
    #assume a priori E(y|x) = 0 for any x
    mean.ml <- crosscov.ml %*% precmat.mlY
    
    
    #sigma^2
    x5.out[i+1] <-rpost_sigma2(mean.ml=mean.ml, z=z)
  }
  out <- cbind(x1.out, x2.out, x3.out, x4.out, x5.out)
  out
}


set.seed(142)
#init_vals<- c(runif(1,.02,.1),runif(1,.02,.4),invgamma::rinvgamma(1,.2,.2))
#init_vals<- c(.02+ (.1-.02)/2, .02+ (.4-.02)/2, 1, 0, .2/(.2+1)) #use mode of IG dist

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")
init_vals<-  as.numeric(res[nrow(res),]); rm(res)

set.seed(625)
pt<-proc.time() # Start Time
#res <- mh.alg(init = init_vals, n.sample = 1e5)
res <- mh.alg(init = init_vals, n.sample = 2e5)
ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3] # End Time to be used in Effective Samples per Second Calculation

#save(res,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")
#save(ptFinal,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/time_homSML.RData")

#save(res,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML2.RData")
#save(ptFinal,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/time_homSML2.RData")

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")
#res1<- res
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML2.RData")
#res2<- res

#res<- rbind(res1,res2)

#nPars=5
#for(i in 1:nPars){
#  plot(density(res[,i]),main=paste("HomML Cheap variable ",i,sep=""),col="black")
#  lines(density(res[1:2e5,i]),col="red")
#}


save(res,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")

#EDIT FROM HERE
################################################################################
################################################################################
####################Performance evaluation######################################
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")

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
plot(1:length(res[,4]),res[,4],xlab="Iteration",ylab="Riverbed elevation error")
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
fit.x.v2 <- scale(fit.x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 

#calculated quantities after fitting homSML emulator
fit.H <- cbind(1, fit.x.test)
fit.crosscov.ml <- ml.crosscov(fit.x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, fit.H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
fit.mean.ml <- fit.crosscov.ml %*% precmat.mlY


mc.mean.ml<- rep(NA,nrow(res))
for(i in 1:nrow(res)){
  
  mc.x.test<- matrix(res[i,],nrow=1,ncol= length(res[i,]))
  
  ##scale validation data accordingly
  mc.x.v2 <- scale(mc.x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
  
  #calculated quantities after fiting homSML emulator
  mc.H <- cbind(1, mc.x.test)
  mc.crosscov.ml <- ml.crosscov(mc.x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, mc.H, var.bc, var.be)
  #assume a priori E(y|x) = 0 for any x
  mc.mean.ml[i] <- mc.crosscov.ml %*% precmat.mlY
}

hist(mc.mean.ml,breaks=20,main="Histogram of predicted Euclidean distance at each MC sample",xlab="Predicted Euclidean Distance")
xtick<-seq(50, 1400, by=50)
axis(side=1, at=xtick, labels = FALSE)

mean(mc.mean.ml); median(mc.mean.ml)

################################################################################
#distribution of predicted euclidean distances from priors of parameters
load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/allParVals.RData")
parVals50m<- as.data.frame(parVals)


prior.mean.ml<- rep(NA,nrow(priorVals))
for(i in 1:nrow(priorVals)){
  
  prior.x.test<- matrix(priorVals[i,],nrow=1,ncol= length(priorVals[i,]))
  
  ##scale validation data accordingly
  prior.x.v2 <- scale(prior.x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
  
  #calculated quantities after fiting homSML emulator
  prior.H <- cbind(1, prior.x.test)
  prior.crosscov.ml <- ml.crosscov(prior.x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, prior.H, var.bc, var.be)
  #assume a priori E(y|x) = 0 for any x
  prior.mean.ml[i] <- prior.crosscov.ml %*% precmat.mlY
}

hist(prior.mean.ml,breaks=20,main="Histogram of predicted Euclidean distance at each prior sample",xlab="Predicted Euclidean Distance")
xtick<-seq(50, 1400, by=50)
axis(side=1, at=xtick, labels = FALSE)

mean(prior.mean.ml); median(prior.mean.ml)


################################################################################
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")
res_HT<- res
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML_HT.RData")

plot(density(res_HT[,1]),main="Channel roughness density",xlab="Channel roughness",col="blue")
lines(density(res[,1]),col="purple")


plot(density(res_HT[,2]),main="Floodplain roughness density",xlab="Floodplain roughness",col="blue")
lines(density(res[,2]),col="purple")


plot(density(res_HT[,3]),main="River width error density",xlab="River width error",col="blue")
lines(density(res[,3]),col="purple")



plot(density(res_HT[,4]),main="Riverbed elevation error density",xlab="Riverbed elevation error",col="blue")
lines(density(res[,4]),col="purple")

