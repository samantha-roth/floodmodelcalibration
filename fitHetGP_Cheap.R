#fit HetGP using only lowres
rm(list=ls())
#setwd(".../FloodingModelCalibrationProject/multires/code")
library(rstan)
library(boot)
library(splines)
library("R.matlab")
source(".../FloodingModelCalibrationProject/multires/code/GPfunctionsOptim.R") #edited code from Kennedy 2020 et al: https://github.com/jcken95/sml-athena
source(".../FloodingModelCalibrationProject/multires/code/hetGPfunctions.R") #edited code from Kennedy 2020 et al: https://github.com/jcken95/sml-athena
hetGP <- stan_model(".../FloodingModelCalibrationProject/multires/code/hetGP.stan") #original code from Kennedy 2020 et al: https://github.com/jcken95/sml-athena

#load parameters
load(".../FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/allParVals.RData")
parVals10m<- as.data.frame(parVals)
load(".../FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/allParVals.RData")
parVals50m<- as.data.frame(parVals)


#load predictions
load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10m.RData")
#50m disaggregated
load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10mfrom50m.RData")
#10m aggregated
#load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics50mfrom10m.RData")


#load true values
#true parameter values
parsTrue<-read.csv(".../FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/parVals/RunTrue_1.csv")


nPars=4

#x data
x.true<- parsTrue[,-(1:2)]

#randomly hold out data to test

################################################################################
##for when testing data is used
set.seed(35)
test.inds<- sample(1:200, 10)

################################################################################
#Random sample of parameters to test on
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

################################################################################
#x.test<- rbind(x.true,parVals10m[test.inds,-1])
#x.test<- rbind(x.true,parVals50m[test.inds,-1])
x.test<- parVals50m[test.inds,-1]
x.hetgp<- parVals50m[-test.inds,-1]
x.c <- parVals50m[-test.inds,-1]
x.e <- parVals10m[-test.inds,-1]

# load data
#y data
#y.test1 <- c(0,metrics.10m$EuclideanDists[test.inds]) #true Euclidean distance between output and reality
#y.test1 <- c(0,metrics.10mfrom50m$EuclideanDists[test.inds]) #true Euclidean distance between output and reality
y.test1 <- metrics.10mfrom50m$EuclideanDists[test.inds] #true Euclidean distance between output and reality
y.hetgp1<- metrics.10mfrom50m$EuclideanDists[-test.inds]
y.c1<- metrics.10mfrom50m$EuclideanDists[-test.inds]
y.e1<- metrics.10m$EuclideanDists[-test.inds]

################################################################################
##no heldout data

#x.test<- as.matrix(x.true)
#x.hetgp<- parVals50m[,-1]
#x.c <- parVals50m[,-1]
#x.e <- parVals10m[,-1]

# load data
#y data
#y.test1 <- 0 #true Euclidean distance between output and reality

#50m disaggregated
#y.hetgp1<- metrics.10mfrom50m$EuclideanDists
#y.c1<- metrics.10mfrom50m$EuclideanDists
#10m aggregated
#y.hetgp1<- metrics.50mfrom10m$EuclideanDists
#y.c1<- metrics.50mfrom10m$EuclideanDists
#y.e1<- metrics.10m$EuclideanDists

################################################################################

##illustration
#################################################################################
#x.test<- x.hetgp
#y.test1<- y.hetgp1
#################################################################################

y.test <- y.test1
y.hetgp <- y.hetgp1
y.c <- y.c1
y.e <- y.e1

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

## standardise inputs (x-mu)/sig

##training data
x.e.std <- scale(x.e)
x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
x.std <- scale(x.hetgp)

##scale validation data accordingly
x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
#x.v1 and x.v2 are the same in my case
m.h.hetgp<- cbind(1, x.hetgp)
m.h.e<- cbind(1, x.e)
m.h.c<- cbind(1, x.c)

#m.h.test<- matrix(c(1,fp.10m.test,ch.10m.test),nrow=1,ncol=length(c(1,fp.10m.test,bspline.ch.10m.test))) #for when only 1 held out obs
m.h.test<- cbind(1,x.test)#for when multiple held out obs

v.h.hetgp<- cbind(1, x.std)

#### fit emulators! ####

### first fit HetGP emulator ####
var.beta <- diag(1, ncol(x.test)+1)
mean.beta <- rep(0, ncol(x.test)+1)

data.hetGP <- list(
  m_p = ncol(x.test)+1, 
  v_p = ncol(x.test)+1,
  N = length(y.hetgp), 
  K = ncol(x.test),
  x = x.std, 
  m_H = m.h.hetgp, 
  v_H = v.h.hetgp,
  y = as.vector(y.hetgp),
  a = rep(1,length(y.hetgp)),
  ## prior
  m_beta_m = rep(0,ncol(x.test)+1), 
  m_beta_s = var.beta,
  m_a_theta = rep(2,ncol(x.test)), 
  m_b_theta = rep(1,ncol(x.test)),
  m_a_sigma = 2, 
  m_b_sigma = 2,
  m_nugget = 0,
  v_beta_m = rep(0, ncol(x.test)+1), 
  v_beta_s = var.beta,
  v_a_theta = rep(2,ncol(x.test)), 
  v_b_theta = rep(1,ncol(x.test)),
  v_a_sigma = 2, 
  v_b_sigma = 2,
  v_nugget_a = 2, 
  v_nugget_b = 2
)
temp <- list()

find.mode <- function(x){
  rstan::optimizing(hetGP, data = data.hetGP, verbose = F, as_vector = F)
}
st.hetgp <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 1)
en.hetgp <- Sys.time() - st.hetgp
en.hetgp
beepr::beep()
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))

het.fit <- temp[[best.emulator]]
pars <- het.fit$par

#save(pars,file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/pars.hetGP.cheap.RData")
#load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/pars.hetGP.cheap.RData")

#save(pars,file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/pars.hetGP.cheap.RData")
#load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/pars.hetGP.cheap.RData")

#save(pars,file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/pars.hetGP.cheap.RData")
load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/pars.hetGP.cheap.RData")



## predict response for hetgp
H <- cbind(1, as.matrix(x.hetgp))

H.new<- cbind(1, as.matrix(x.test))

covmat.hetgp <- H %*% var.beta %*% t(H) + 
  cov.x1.x2(x.std, x.std, pars$m_sigma, pars$m_theta, 2) + 
  diag(exp(pars$logLambda))

precmat.hetgp <- chol2inv(chol(covmat.hetgp))

crosscov <- gp.crosscov(x.v1, x.std, pars$m_sigma, pars$m_theta, H.new, H, var.beta )

mean.hetgp <- crosscov %*% (precmat.hetgp %*% (y.hetgp - H %*% mean.beta))

var.hetgp <- H.new %*% var.beta %*% t(H.new) + cov.x1.x2(x.v1, x.v1, pars$m_sigma, pars$m_theta, 2)
var.hetgp <- var.hetgp - crosscov %*% precmat.hetgp %*% t(crosscov)

H.l <- cbind(1, x.std)
H.l.new <- cbind(1, x.v1)

covmat.hetgp.lambda <- H.l %*% var.beta %*% t(H.l) + 
  cov.x1.x2(x.std, x.std, pars$v_sigma, pars$v_theta, 2) + 
  diag(pars$v_nugget, dim(x.std)[1])

precmat.hetgp.lambda <- chol2inv(chol(covmat.hetgp.lambda))

crosscov.hetgp.lambda <- gp.crosscov(x.v1, x.std, pars$v_sigma, pars$v_theta,H.l.new, H.l, var.beta)
lambda.hetgp <- exp( crosscov.hetgp.lambda %*% precmat.hetgp.lambda %*%(pars$logLambda - as.vector(H.l%*%mean.beta) ))
var.hetgp.full <- var.hetgp + diag(lambda.hetgp)

plot(mean.hetgp, y.test)

plot(x.test[,1],mean.hetgp) 
plot(x.test[,1],lambda.hetgp)

plot(x.test[,2],mean.hetgp) 
plot(x.test[,2],lambda.hetgp)

plot(x.test[,3],mean.hetgp) 
plot(x.test[,3],lambda.hetgp)

plot(x.test[,4],mean.hetgp) 
plot(x.test[,4],lambda.hetgp)


MSE.hetgpcheap <- MSE(y.test, mean.hetgp)
mean.hetgpcheap<- mean.hetgp
lambda.hetgpcheap<- lambda.hetgp

x.test.cheap<- x.test
y.test.cheap<- y.test

save(mean.hetgpcheap, lambda.hetgp, x.test.cheap, y.test.cheap, MSE.hetgpcheap, file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/hetGPCheap_calculated_quantities.RData")


save(mean.hetgp, lambda.hetgp, MSE.hetgp, file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/hetGPCheap_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/hetGPCheap_calculated_quantities.RData")
