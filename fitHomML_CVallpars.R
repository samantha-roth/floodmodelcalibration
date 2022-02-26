## this is the one to use
rm(list=ls())
setwd("C:/FloodingModelCalibrationProject/sml-athena-main")
library(rstan)
library(boot)
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
#randomly hold out data to test

################################################################################
#for when 5 sets of samples of size 20 are used

#for when testing data is used
set.seed(14)
randOrderInds<- sample(1:200, 200)
groups<- rep(1:10, each=20)

################################################################################

#first test set
test.inds<- randOrderInds[1:20]

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

################################################################################
############################ fit SML emulator ##################################
################################################################################

## need to get the data into shape!
n.ml <- length(y.e)

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2); n.e <- length(y.e2)
x.c.std.rev <- x.c.std[n.c:1,]
x.e.std.rev <- x.e.std[n.ml:1,]

#pairs(cbind(y.c2, x.c.std))
#pairs(cbind(y.e2, x.e.std))
## reverse everything to get the order correct ...
n.c <- length(y.c2); n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

x.c.rev<- as.matrix(x.c[n.c:1,])
x.e.rev<- as.matrix(x.e[n.e:1,])

m.h <- cbind(cbind(1, x.c.rev), rbind(matrix(0, ncol = (ncol(x.test)+1), nrow = n.c - n.e) ,cbind(1, x.e.rev)))
#tail(m.h)
var.beta <- diag(1, ncol(m.h))

ml.data <- list(
  
  ## data ##
  
  m_p = (ncol(x.test)+1),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  
  m_beta_m = rep(0, (ncol(x.test)+1)), 
  m_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget_a = 2, m_nugget_b = 2,
  
  c_beta_m = rep(0, (ncol(x.test)+1)), 
  c_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  c_a_theta = rep(2, ncol(x.test)), 
  c_b_theta = rep(1, ncol(x.test)),
  c_a_sigma = 2, 
  c_b_sigma = 2,
  c_nugget_a = 2, 
  c_nugget_b = 2,
  
  m_rho = 1, 
  s_rho = 1/3
  
)


temp <- list()

## fit multilevel GP

find.mode <- function(x){
  rstan::optimizing(homSML, data = ml.data, verbose = F, as_vector = F)
}
st1 <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 1)
en1 <- Sys.time()
en1-st1
beepr::beep(4)
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))
c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value)

ml.fit <-  temp[[best.emulator]]
pars <- ml.fit$par


save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV1pars.homML.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV1pars.homML.RData")


designmat <- matrix(0, ncol=ncol(m.h), nrow=length(c(y.c, y.e)))
designmat[1:length(y.c),1:(ncol(x.test)+1)] <- cbind(1, x.c.rev)
designmat[-(1:length(y.c)), (1:(ncol(x.test)+1))] <- pars$rho*cbind(1, x.e.rev)
designmat[-(1:length(y.c)), -(1:(ncol(x.test)+1))] <- cbind(1, x.e.rev)
dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
precmat.ml <- chol2inv(chol(dat.covmat))
rm(dat.covmat)
H.c <- cbind(1, x.c.rev)
H.e <- cbind(1, as.matrix(x.c[n.e:1,]))
H <- cbind(1, as.matrix(x.test))

var.bc <- var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)]
var.be <- var.beta[(ncol(x.test)+2):(ncol(m.h)), (ncol(x.test)+2):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
var.ml.full <- var.ml + diag(rep(pars$m_nugget,nrow(var.ml)))

#### compare the fits ####

## score & MSE

MSE.hommlCV1 <- MSE(y.test, mean.ml)

y.test.hommlCV1<- y.test
x.test.hommlCV1<- x.test
mean.hommlCV1<- mean.ml
#save values of metrics

save(MSE.hommlCV1, y.test.hommlCV1, x.test.hommlCV1, mean.hommlCV1,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV1homML_calculated_quantities.RData")

################################################################################
################################################################################

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

################################################################################
############################ fit SML emulator ##################################
################################################################################

## need to get the data into shape!
n.ml <- length(y.e)

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2); n.e <- length(y.e2)
x.c.std.rev <- x.c.std[n.c:1,]
x.e.std.rev <- x.e.std[n.ml:1,]

#pairs(cbind(y.c2, x.c.std))
#pairs(cbind(y.e2, x.e.std))
## reverse everything to get the order correct ...
n.c <- length(y.c2); n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

x.c.rev<- as.matrix(x.c[n.c:1,])
x.e.rev<- as.matrix(x.e[n.e:1,])

m.h <- cbind(cbind(1, x.c.rev), rbind(matrix(0, ncol = (ncol(x.test)+1), nrow = n.c - n.e) ,cbind(1, x.e.rev)))
#tail(m.h)
var.beta <- diag(1, ncol(m.h))

ml.data <- list(
  
  ## data ##
  
  m_p = (ncol(x.test)+1),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  
  m_beta_m = rep(0, (ncol(x.test)+1)), 
  m_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget_a = 2, m_nugget_b = 2,
  
  c_beta_m = rep(0, (ncol(x.test)+1)), 
  c_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  c_a_theta = rep(2, ncol(x.test)), 
  c_b_theta = rep(1, ncol(x.test)),
  c_a_sigma = 2, 
  c_b_sigma = 2,
  c_nugget_a = 2, 
  c_nugget_b = 2,
  
  m_rho = 1, 
  s_rho = 1/3
  
)


temp <- list()

## fit multilevel GP

find.mode <- function(x){
  rstan::optimizing(homSML, data = ml.data, verbose = F, as_vector = F)
}
st1 <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 1)
en1 <- Sys.time()
en1-st1
beepr::beep(4)
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))
c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value)

ml.fit <-  temp[[best.emulator]]
pars <- ml.fit$par


save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2pars.homML.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2pars.homML.RData")


designmat <- matrix(0, ncol=ncol(m.h), nrow=length(c(y.c, y.e)))
designmat[1:length(y.c),1:(ncol(x.test)+1)] <- cbind(1, x.c.rev)
designmat[-(1:length(y.c)), (1:(ncol(x.test)+1))] <- pars$rho*cbind(1, x.e.rev)
designmat[-(1:length(y.c)), -(1:(ncol(x.test)+1))] <- cbind(1, x.e.rev)
dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
precmat.ml <- chol2inv(chol(dat.covmat))
rm(dat.covmat)
H.c <- cbind(1, x.c.rev)
H.e <- cbind(1, as.matrix(x.c[n.e:1,]))
H <- cbind(1, as.matrix(x.test))

var.bc <- var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)]
var.be <- var.beta[(ncol(x.test)+2):(ncol(m.h)), (ncol(x.test)+2):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
var.ml.full <- var.ml + diag(rep(pars$m_nugget,nrow(var.ml)))

#### compare the fits ####

## score & MSE

MSE.hommlCV2 <- MSE(y.test, mean.ml)

y.test.hommlCV2<- y.test
x.test.hommlCV2<- x.test
mean.hommlCV2<- mean.ml
#save values of metrics

save(MSE.hommlCV2, y.test.hommlCV2, x.test.hommlCV2, mean.hommlCV2,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV2homML_calculated_quantities.RData")

################################################################################
################################################################################
#third test set
test.inds<- randOrderInds[41:60]

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

################################################################################
############################ fit SML emulator ##################################
################################################################################

## need to get the data into shape!
n.ml <- length(y.e)

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2); n.e <- length(y.e2)
x.c.std.rev <- x.c.std[n.c:1,]
x.e.std.rev <- x.e.std[n.ml:1,]

#pairs(cbind(y.c2, x.c.std))
#pairs(cbind(y.e2, x.e.std))
## reverse everything to get the order correct ...
n.c <- length(y.c2); n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

x.c.rev<- as.matrix(x.c[n.c:1,])
x.e.rev<- as.matrix(x.e[n.e:1,])

m.h <- cbind(cbind(1, x.c.rev), rbind(matrix(0, ncol = (ncol(x.test)+1), nrow = n.c - n.e) ,cbind(1, x.e.rev)))
#tail(m.h)
var.beta <- diag(1, ncol(m.h))

ml.data <- list(
  
  ## data ##
  
  m_p = (ncol(x.test)+1),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  
  m_beta_m = rep(0, (ncol(x.test)+1)), 
  m_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget_a = 2, m_nugget_b = 2,
  
  c_beta_m = rep(0, (ncol(x.test)+1)), 
  c_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  c_a_theta = rep(2, ncol(x.test)), 
  c_b_theta = rep(1, ncol(x.test)),
  c_a_sigma = 2, 
  c_b_sigma = 2,
  c_nugget_a = 2, 
  c_nugget_b = 2,
  
  m_rho = 1, 
  s_rho = 1/3
  
)


temp <- list()

## fit multilevel GP

find.mode <- function(x){
  rstan::optimizing(homSML, data = ml.data, verbose = F, as_vector = F)
}
st1 <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 1)
en1 <- Sys.time()
en1-st1
beepr::beep(4)
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))
c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value)

ml.fit <-  temp[[best.emulator]]
pars <- ml.fit$par


save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV3pars.homML.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV3pars.homML.RData")


designmat <- matrix(0, ncol=ncol(m.h), nrow=length(c(y.c, y.e)))
designmat[1:length(y.c),1:(ncol(x.test)+1)] <- cbind(1, x.c.rev)
designmat[-(1:length(y.c)), (1:(ncol(x.test)+1))] <- pars$rho*cbind(1, x.e.rev)
designmat[-(1:length(y.c)), -(1:(ncol(x.test)+1))] <- cbind(1, x.e.rev)
dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
precmat.ml <- chol2inv(chol(dat.covmat))
rm(dat.covmat)
H.c <- cbind(1, x.c.rev)
H.e <- cbind(1, as.matrix(x.c[n.e:1,]))
H <- cbind(1, as.matrix(x.test))

var.bc <- var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)]
var.be <- var.beta[(ncol(x.test)+2):(ncol(m.h)), (ncol(x.test)+2):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
var.ml.full <- var.ml + diag(rep(pars$m_nugget,nrow(var.ml)))

#### compare the fits ####

## score & MSE

MSE.hommlCV3 <- MSE(y.test, mean.ml)

y.test.hommlCV3<- y.test
x.test.hommlCV3<- x.test
mean.hommlCV3<- mean.ml
#save values of metrics

save(MSE.hommlCV3, y.test.hommlCV3, x.test.hommlCV3, mean.hommlCV3,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV3homML_calculated_quantities.RData")

################################################################################
################################################################################

#4th test set
test.inds<- randOrderInds[61:80]

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

################################################################################
############################ fit SML emulator ##################################
################################################################################

## need to get the data into shape!
n.ml <- length(y.e)

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2); n.e <- length(y.e2)
x.c.std.rev <- x.c.std[n.c:1,]
x.e.std.rev <- x.e.std[n.ml:1,]

#pairs(cbind(y.c2, x.c.std))
#pairs(cbind(y.e2, x.e.std))
## reverse everything to get the order correct ...
n.c <- length(y.c2); n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

x.c.rev<- as.matrix(x.c[n.c:1,])
x.e.rev<- as.matrix(x.e[n.e:1,])

m.h <- cbind(cbind(1, x.c.rev), rbind(matrix(0, ncol = (ncol(x.test)+1), nrow = n.c - n.e) ,cbind(1, x.e.rev)))
#tail(m.h)
var.beta <- diag(1, ncol(m.h))

ml.data <- list(
  
  ## data ##
  
  m_p = (ncol(x.test)+1),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  
  m_beta_m = rep(0, (ncol(x.test)+1)), 
  m_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget_a = 2, m_nugget_b = 2,
  
  c_beta_m = rep(0, (ncol(x.test)+1)), 
  c_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  c_a_theta = rep(2, ncol(x.test)), 
  c_b_theta = rep(1, ncol(x.test)),
  c_a_sigma = 2, 
  c_b_sigma = 2,
  c_nugget_a = 2, 
  c_nugget_b = 2,
  
  m_rho = 1, 
  s_rho = 1/3
  
)


temp <- list()

## fit multilevel GP

find.mode <- function(x){
  rstan::optimizing(homSML, data = ml.data, verbose = F, as_vector = F)
}
st1 <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 1)
en1 <- Sys.time()
en1-st1
beepr::beep(4)
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))
c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value)

ml.fit <-  temp[[best.emulator]]
pars <- ml.fit$par


save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV4pars.homML.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV4pars.homML.RData")


designmat <- matrix(0, ncol=ncol(m.h), nrow=length(c(y.c, y.e)))
designmat[1:length(y.c),1:(ncol(x.test)+1)] <- cbind(1, x.c.rev)
designmat[-(1:length(y.c)), (1:(ncol(x.test)+1))] <- pars$rho*cbind(1, x.e.rev)
designmat[-(1:length(y.c)), -(1:(ncol(x.test)+1))] <- cbind(1, x.e.rev)
dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
precmat.ml <- chol2inv(chol(dat.covmat))
rm(dat.covmat)
H.c <- cbind(1, x.c.rev)
H.e <- cbind(1, as.matrix(x.c[n.e:1,]))
H <- cbind(1, as.matrix(x.test))

var.bc <- var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)]
var.be <- var.beta[(ncol(x.test)+2):(ncol(m.h)), (ncol(x.test)+2):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
var.ml.full <- var.ml + diag(rep(pars$m_nugget,nrow(var.ml)))

#### compare the fits ####

## score & MSE

MSE.hommlCV4 <- MSE(y.test, mean.ml)

y.test.hommlCV4<- y.test
x.test.hommlCV4<- x.test
mean.hommlCV4<- mean.ml
#save values of metrics

save(MSE.hommlCV4, y.test.hommlCV4, x.test.hommlCV4, mean.hommlCV4,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV4homML_calculated_quantities.RData")

################################################################################
################################################################################
#5th test set
test.inds<- randOrderInds[81:100]

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

################################################################################
############################ fit SML emulator ##################################
################################################################################

## need to get the data into shape!
n.ml <- length(y.e)

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2); n.e <- length(y.e2)
x.c.std.rev <- x.c.std[n.c:1,]
x.e.std.rev <- x.e.std[n.ml:1,]

#pairs(cbind(y.c2, x.c.std))
#pairs(cbind(y.e2, x.e.std))
## reverse everything to get the order correct ...
n.c <- length(y.c2); n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

x.c.rev<- as.matrix(x.c[n.c:1,])
x.e.rev<- as.matrix(x.e[n.e:1,])

m.h <- cbind(cbind(1, x.c.rev), rbind(matrix(0, ncol = (ncol(x.test)+1), nrow = n.c - n.e) ,cbind(1, x.e.rev)))
#tail(m.h)
var.beta <- diag(1, ncol(m.h))

ml.data <- list(
  
  ## data ##
  
  m_p = (ncol(x.test)+1),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  
  m_beta_m = rep(0, (ncol(x.test)+1)), 
  m_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget_a = 2, m_nugget_b = 2,
  
  c_beta_m = rep(0, (ncol(x.test)+1)), 
  c_beta_s = var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)],
  c_a_theta = rep(2, ncol(x.test)), 
  c_b_theta = rep(1, ncol(x.test)),
  c_a_sigma = 2, 
  c_b_sigma = 2,
  c_nugget_a = 2, 
  c_nugget_b = 2,
  
  m_rho = 1, 
  s_rho = 1/3
  
)


temp <- list()

## fit multilevel GP

find.mode <- function(x){
  rstan::optimizing(homSML, data = ml.data, verbose = F, as_vector = F)
}
st1 <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 1)
en1 <- Sys.time()
en1-st1
beepr::beep(4)
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))
c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value)

ml.fit <-  temp[[best.emulator]]
pars <- ml.fit$par


save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV5pars.homML.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV5pars.homML.RData")

designmat <- matrix(0, ncol=ncol(m.h), nrow=length(c(y.c, y.e)))
designmat[1:length(y.c),1:(ncol(x.test)+1)] <- cbind(1, x.c.rev)
designmat[-(1:length(y.c)), (1:(ncol(x.test)+1))] <- pars$rho*cbind(1, x.e.rev)
designmat[-(1:length(y.c)), -(1:(ncol(x.test)+1))] <- cbind(1, x.e.rev)
dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
precmat.ml <- chol2inv(chol(dat.covmat))
rm(dat.covmat)
H.c <- cbind(1, x.c.rev)
H.e <- cbind(1, as.matrix(x.c[n.e:1,]))
H <- cbind(1, as.matrix(x.test))

var.bc <- var.beta[1:(ncol(x.test)+1),1:(ncol(x.test)+1)]
var.be <- var.beta[(ncol(x.test)+2):(ncol(m.h)), (ncol(x.test)+2):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
var.ml.full <- var.ml + diag(rep(pars$m_nugget,nrow(var.ml)))

#### compare the fits ####

## score & MSE

MSE.hommlCV5 <- MSE(y.test, mean.ml)

y.test.hommlCV5<- y.test
x.test.hommlCV5<- x.test
mean.hommlCV5<- mean.ml
#save values of metrics

save(MSE.hommlCV5, y.test.hommlCV5, x.test.hommlCV5, mean.hommlCV5,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV5homML_calculated_quantities.RData")

################################################################################
################################################################################