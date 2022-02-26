## this is the one to use
rm(list=ls())
setwd("C:/FloodingModelCalibrationProject/sml-athena-main")
library(raster)
library(rstan)
#library(boot)
library("R.matlab")
source("GPfunctionsOptim.R")
source("hetGPfunctions.R")
#hetGP <- stan_model("hetGP.stan")
homGP<- stan_model("GP.stan")
#sml <- stan_model("het-SML.stan")
homSML<- stan_model("hom-SML.stan")

#load parameters
load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/allParVals.RData")
parVals10m<- as.data.frame(parVals)
load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/allParVals.RData")
parVals50m<- as.data.frame(parVals)

#load predictions
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10m.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10mfrom50m.RData")

#load true values
#true parameter values
parsTrue<- read.csv("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/parVals/RunTrue_1.csv")


nPars=4

#x data
x.true<- parsTrue[,-(1:2)]

#randomly hold out data to test

################################################################################
#for when 5 sets of samples of size 20 are used

#for when testing data is used
set.seed(14)
randOrderInds<- sample(1:200, 200)
groups<- rep(1:10, each=20)

#first test set
test.inds<- randOrderInds[21:40]

################################################################################
#x.test<- rbind(x.true,parVals10m[test.inds,-1])
x.test<- parVals10m[test.inds,-1]
x.gp<- parVals10m[-test.inds,-1]
x.c <- parVals50m[-test.inds,-1]
x.e <- parVals10m[-test.inds,-1]

# load data
#y data
#y.test1 <- c(0,metrics.10m$EuclideanDists[test.inds]) #true Euclidean distance between output and reality
y.test1 <- metrics.10m$EuclideanDists[test.inds]
y.gp1<- metrics.10m$EuclideanDists[-test.inds]
y.c1<- metrics.10mfrom50m$EuclideanDists[-test.inds]
y.e1<- metrics.10m$EuclideanDists[-test.inds]

################################################################################

y.test <- y.test1
y.gp <- y.gp1
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
x.std <- scale(x.gp)

##scale validation data accordingly
x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))
#x.v1 and x.v2 are the same in my case

#### fit emulators! ####

### first fit HomGP emulator ####
var.beta <- diag(1, ncol(x.test)+1)
mean.beta <- rep(0, ncol(x.test)+1)
data.GP <- list(
  m_p = ncol(x.test)+1, 
  N = length(y.gp), 
  K = ncol(x.test),
  x = x.std, 
  m_H = cbind(1, x.gp), 
  y = as.vector(y.gp),
  a = rep(1,length(y.gp)),
  ## prior
  m_beta_m = rep(0,ncol(x.test)+1), 
  m_beta_s = var.beta,
  m_a_theta = rep(2,ncol(x.test)), 
  m_b_theta = rep(1,ncol(x.test)),
  m_a_sigma = 2, 
  m_b_sigma = 2,
  m_nugget_a = 2, 
  m_nugget_b = 2
)
temp <- list()

find.mode <- function(x){
  rstan::optimizing(homGP, data = data.GP, verbose = F, as_vector = F)
}
st.gp <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 1)
en.gp <- Sys.time() - st.gp
en.gp
beepr::beep()
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))

gp.fit <- temp[[best.emulator]]


pars <- gp.fit$par

save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2pars.homGP.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2pars.homGP.RData")


### predict response for gp
H <- cbind(1, as.matrix(x.gp))

H.new <- cbind(1, as.matrix(x.test))

covmat.gp <- H %*% var.beta %*% t(H) + 
  cov.x1.x2(x.std, x.std, pars$m_sigma, pars$m_theta, 2) + diag(rep(pars$m_nugget),length(y.gp))

precmat.gp <- chol2inv(chol(covmat.gp))

crosscov <- gp.crosscov(x.v1, x.std, pars$m_sigma, pars$m_theta, H.new, H, var.beta )

mean.gp <- crosscov %*% (precmat.gp %*% (y.gp - H %*% mean.beta))

var.gp <- H.new %*% var.beta %*% t(H.new) + cov.x1.x2(x.v1, x.v1, pars$m_sigma, pars$m_theta, 2)
var.gp <- var.gp - crosscov %*% precmat.gp %*% t(crosscov)

plot(mean.gp, y.test)
mean((mean.gp - y.test)^2)

MSE.gpCV2<- MSE(y.test,mean.gp)

y.test.gpCV2<- y.test
x.test.gpCV2<- x.test
mean.gpCV2<- mean.gp

save(mean.gpCV2, x.test.gpCV2, y.test.gpCV2, MSE.gpCV2,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2homGP_calculated_quantities.RData")
