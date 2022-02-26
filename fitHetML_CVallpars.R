#fit HetGP and HetML to 4 parameter sample
rm(list=ls())
#setwd(".../FloodingModelCalibrationProject/multires/code")
library(rstan)
library(boot)
library(splines)
library("R.matlab")
source(".../FloodingModelCalibrationProject/multires/code/GPfunctionsOptim.R") #edited code from Kennedy et al 2020: https://github.com/jcken95/sml-athena
source(".../FloodingModelCalibrationProject/multires/code/hetGPfunctions.R") #edited code from Kennedy et al 2020: https://github.com/jcken95/sml-athena
#hetGP <- stan_model(".../FloodingModelCalibrationProject/multires/code/hetGP.stan") #original code from Kennedy et al 2020: https://github.com/jcken95/sml-athena
sml <- stan_model(".../FloodingModelCalibrationProject/multires/code/het-SML.stan") #original code from Kennedy et al 2020: https://github.com/jcken95/sml-athena

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
#for when 5 sets of samples of size 20 are used

#for when testing data is used
set.seed(14)
randOrderInds<- sample(1:200, 200)
groups<- rep(1:10, each=20)


################################################################################
################################################################################

#first test set
test.inds<- randOrderInds[1:20]

################################################################################
#x.test<- rbind(x.true,parVals10m[test.inds,-1])
x.test<- rbind(parVals10m[test.inds,-1])
x.hetgp<- parVals10m[-test.inds,-1]
x.c <- parVals50m[-test.inds,-1]
x.e <- parVals10m[-test.inds,-1]

# load data
#y data
#y.test1 <- c(0,metrics.10m$EuclideanDists[test.inds]) #true Euclidean distance between output and reality
y.test1 <- metrics.10m$EuclideanDists[test.inds]
y.hetgp1<- metrics.10m$EuclideanDists[-test.inds]
y.c1<- metrics.10mfrom50m$EuclideanDists[-test.inds]
y.e1<- metrics.10m$EuclideanDists[-test.inds]
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
n.c <- length(y.c2)
n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

#create H
m.h.c.rev<- as.matrix(m.h.c[n.c:1,])

m.h.0<- matrix(0, ncol = ncol(m.h.c.rev), nrow = nrow(m.h.c.rev) - n.e)
m.h.e.rev<- as.matrix(m.h.e[n.e:1,])
m.h.2<- rbind(m.h.0, m.h.e.rev)

m.h <- cbind(m.h.c.rev, m.h.2)

v.h.rev<- cbind(1, x.e.std.rev)

#tail(m.h)
var.beta <- diag(1, ncol(m.h))
var.beta.m.h.c<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]
var.beta.m.h.v<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]

ml.data <- list(
  ## data ##
  m_p = ncol(m.h.e), 
  v_p = ncol(v.h.rev),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h, 
  v_H = v.h.rev,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  m_beta_m = rep(0, ncol(m.h.e)), 
  m_beta_s = var.beta.m.h.c,
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, 
  m_b_sigma = 2,
  m_nugget = 0,
  
  v_beta_m = rep(0, ncol(m.h.e)),
  v_beta_s = var.beta.m.h.v,
  v_a_theta = rep(2, ncol(x.test)), 
  v_b_theta = rep(1, ncol(x.test)),
  v_a_sigma = 2, 
  v_b_sigma = 2,
  v_nugget_a = 2, 
  v_nugget_b = 2,
  
  c_beta_m = rep(0, ncol(m.h.c)), 
  c_beta_s = var.beta.m.h.c,
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
  rstan::optimizing(sml, data = ml.data, verbose = F, as_vector = F)
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

save(pars,file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV1pars.hetML.RData")

load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV1pars.hetML.RData")

designmat <- matrix(0, ncol=ncol(m.h), nrow=n.c+n.e)
designmat[1:length(y.c),1:(ncol(m.h.c.rev))] <- m.h.c.rev
designmat[-(1:length(y.c)), (1:(ncol(m.h.e.rev)))] <- pars$rho*m.h.e.rev
designmat[-(1:length(y.c)), -(1:(ncol(m.h.e.rev)))] <- m.h.e.rev
#designmat is what is referred to as H on page 11 of the revised SML paper

dat.covmat <- dataCovmat(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
#dat.covmat = Var(Y|all parameters) + H%*%var.beta%*%t(H) on page 10 of revised SML paper
precmat.ml <- chol2inv(chol(dat.covmat))
#precmat.ml = (Var(Y|all parameters) + H%*%var.beta%*%t(H))^-1
rm(dat.covmat)

H.c <- as.matrix(m.h.c.rev)
H.e <- as.matrix(m.h.e.rev)
H.ml.lambda <- as.matrix(v.h.rev)
H <- as.matrix(m.h.test)

var.bc <- var.beta[1:ncol(m.h.c),1:ncol(m.h.c)]
var.be <- var.beta[(ncol(m.h.c)+1):(ncol(m.h)), (ncol(m.h.c)+1):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
covmat.ml.lambda <- H.ml.lambda %*% var.beta[1:ncol(m.h.e),1:ncol(m.h.e)] %*% t(H.ml.lambda) + cov.x1.x2(x.e.std.rev, x.e.std.rev, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.e.std.rev)[1])
precmat.ml.lambda <- chol2inv(chol(covmat.ml.lambda))
crosscov.ml.lambda <- gp.crosscov(x.v2, x.e.std.rev, pars$v_sigma, pars$v_theta, H, H.ml.lambda, var.beta[1:ncol(m.h.e),1:ncol(m.h.e)])
lambda.ml <- exp( crosscov.ml.lambda %*%( precmat.ml.lambda %*% pars$logLambda) )
var.ml.full <- var.ml + diag(lambda.ml)

MSE.mlCV1 <- MSE(y.test, mean.ml)

x.test.hetmlCV1<- x.test
y.test.hetmlCV1<- y.test

mean.mlCV1= mean.ml
lambda.mlCV1= lambda.ml

save(MSE.mlCV1,mean.mlCV1,lambda.mlCV1,x.test.hetmlCV1,y.test.hetmlCV1,
     file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV1hetml_calculated_quantities.RData")


################################################################################
################################################################################

#2nd test set
test.inds<- randOrderInds[21:40]

################################################################################
#x.test<- rbind(x.true,parVals10m[test.inds,-1])
x.test<- rbind(parVals10m[test.inds,-1])
x.hetgp<- parVals10m[-test.inds,-1]
x.c <- parVals50m[-test.inds,-1]
x.e <- parVals10m[-test.inds,-1]

# load data
#y data
#y.test1 <- c(0,metrics.10m$EuclideanDists[test.inds]) #true Euclidean distance between output and reality
y.test1 <- metrics.10m$EuclideanDists[test.inds]
y.hetgp1<- metrics.10m$EuclideanDists[-test.inds]
y.c1<- metrics.10mfrom50m$EuclideanDists[-test.inds]
y.e1<- metrics.10m$EuclideanDists[-test.inds]
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
n.c <- length(y.c2)
n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

#create H
m.h.c.rev<- as.matrix(m.h.c[n.c:1,])

m.h.0<- matrix(0, ncol = ncol(m.h.c.rev), nrow = nrow(m.h.c.rev) - n.e)
m.h.e.rev<- as.matrix(m.h.e[n.e:1,])
m.h.2<- rbind(m.h.0, m.h.e.rev)

m.h <- cbind(m.h.c.rev, m.h.2)

v.h.rev<- cbind(1, x.e.std.rev)

#tail(m.h)
var.beta <- diag(1, ncol(m.h))
var.beta.m.h.c<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]
var.beta.m.h.v<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]

ml.data <- list(
  ## data ##
  m_p = ncol(m.h.e), 
  v_p = ncol(v.h.rev),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h, 
  v_H = v.h.rev,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  m_beta_m = rep(0, ncol(m.h.e)), 
  m_beta_s = var.beta.m.h.c,
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, 
  m_b_sigma = 2,
  m_nugget = 0,
  
  v_beta_m = rep(0, ncol(m.h.e)),
  v_beta_s = var.beta.m.h.v,
  v_a_theta = rep(2, ncol(x.test)), 
  v_b_theta = rep(1, ncol(x.test)),
  v_a_sigma = 2, 
  v_b_sigma = 2,
  v_nugget_a = 2, 
  v_nugget_b = 2,
  
  c_beta_m = rep(0, ncol(m.h.c)), 
  c_beta_s = var.beta.m.h.c,
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
  rstan::optimizing(sml, data = ml.data, verbose = F, as_vector = F)
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

save(pars,file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2pars.hetML.RData")

load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2pars.hetML.RData")

designmat <- matrix(0, ncol=ncol(m.h), nrow=n.c+n.e)
designmat[1:length(y.c),1:(ncol(m.h.c.rev))] <- m.h.c.rev
designmat[-(1:length(y.c)), (1:(ncol(m.h.e.rev)))] <- pars$rho*m.h.e.rev
designmat[-(1:length(y.c)), -(1:(ncol(m.h.e.rev)))] <- m.h.e.rev
#designmat is what is referred to as H on page 11 of the revised SML paper

dat.covmat <- dataCovmat(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
#dat.covmat = Var(Y|all parameters) + H%*%var.beta%*%t(H) on page 10 of revised SML paper
precmat.ml <- chol2inv(chol(dat.covmat))
#precmat.ml = (Var(Y|all parameters) + H%*%var.beta%*%t(H))^-1
rm(dat.covmat)

H.c <- as.matrix(m.h.c.rev)
H.e <- as.matrix(m.h.e.rev)
H.ml.lambda <- as.matrix(v.h.rev)
H <- as.matrix(m.h.test)

var.bc <- var.beta[1:ncol(m.h.c),1:ncol(m.h.c)]
var.be <- var.beta[(ncol(m.h.c)+1):(ncol(m.h)), (ncol(m.h.c)+1):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
covmat.ml.lambda <- H.ml.lambda %*% var.beta[1:ncol(m.h.e),1:ncol(m.h.e)] %*% t(H.ml.lambda) + cov.x1.x2(x.e.std.rev, x.e.std.rev, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.e.std.rev)[1])
precmat.ml.lambda <- chol2inv(chol(covmat.ml.lambda))
crosscov.ml.lambda <- gp.crosscov(x.v2, x.e.std.rev, pars$v_sigma, pars$v_theta, H, H.ml.lambda, var.beta[1:ncol(m.h.e),1:ncol(m.h.e)])
lambda.ml <- exp( crosscov.ml.lambda %*%( precmat.ml.lambda %*% pars$logLambda) )
var.ml.full <- var.ml + diag(lambda.ml)

MSE.mlCV2 <- MSE(y.test, mean.ml)

x.test.hetmlCV2<- x.test
y.test.hetmlCV2<- y.test

mean.mlCV2= mean.ml
lambda.mlCV2= lambda.ml

save(MSE.mlCV2,mean.mlCV2,lambda.mlCV2,x.test.hetmlCV2,y.test.hetmlCV2,
     file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2hetml_calculated_quantities.RData")

################################################################################
################################################################################

#3rd test set
test.inds<- randOrderInds[41:60]

################################################################################
#x.test<- rbind(x.true,parVals10m[test.inds,-1])
x.test<- rbind(parVals10m[test.inds,-1])
x.hetgp<- parVals10m[-test.inds,-1]
x.c <- parVals50m[-test.inds,-1]
x.e <- parVals10m[-test.inds,-1]

# load data
#y data
#y.test1 <- c(0,metrics.10m$EuclideanDists[test.inds]) #true Euclidean distance between output and reality
y.test1 <- metrics.10m$EuclideanDists[test.inds]
y.hetgp1<- metrics.10m$EuclideanDists[-test.inds]
y.c1<- metrics.10mfrom50m$EuclideanDists[-test.inds]
y.e1<- metrics.10m$EuclideanDists[-test.inds]
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
n.c <- length(y.c2)
n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

#create H
m.h.c.rev<- as.matrix(m.h.c[n.c:1,])

m.h.0<- matrix(0, ncol = ncol(m.h.c.rev), nrow = nrow(m.h.c.rev) - n.e)
m.h.e.rev<- as.matrix(m.h.e[n.e:1,])
m.h.2<- rbind(m.h.0, m.h.e.rev)

m.h <- cbind(m.h.c.rev, m.h.2)

v.h.rev<- cbind(1, x.e.std.rev)

#tail(m.h)
var.beta <- diag(1, ncol(m.h))
var.beta.m.h.c<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]
var.beta.m.h.v<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]

ml.data <- list(
  ## data ##
  m_p = ncol(m.h.e), 
  v_p = ncol(v.h.rev),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h, 
  v_H = v.h.rev,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  m_beta_m = rep(0, ncol(m.h.e)), 
  m_beta_s = var.beta.m.h.c,
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, 
  m_b_sigma = 2,
  m_nugget = 0,
  
  v_beta_m = rep(0, ncol(m.h.e)),
  v_beta_s = var.beta.m.h.v,
  v_a_theta = rep(2, ncol(x.test)), 
  v_b_theta = rep(1, ncol(x.test)),
  v_a_sigma = 2, 
  v_b_sigma = 2,
  v_nugget_a = 2, 
  v_nugget_b = 2,
  
  c_beta_m = rep(0, ncol(m.h.c)), 
  c_beta_s = var.beta.m.h.c,
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
  rstan::optimizing(sml, data = ml.data, verbose = F, as_vector = F)
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

save(pars,file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV3pars.hetML.RData")

load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV3pars.hetML.RData")

designmat <- matrix(0, ncol=ncol(m.h), nrow=n.c+n.e)
designmat[1:length(y.c),1:(ncol(m.h.c.rev))] <- m.h.c.rev
designmat[-(1:length(y.c)), (1:(ncol(m.h.e.rev)))] <- pars$rho*m.h.e.rev
designmat[-(1:length(y.c)), -(1:(ncol(m.h.e.rev)))] <- m.h.e.rev
#designmat is what is referred to as H on page 11 of the revised SML paper

dat.covmat <- dataCovmat(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
#dat.covmat = Var(Y|all parameters) + H%*%var.beta%*%t(H) on page 10 of revised SML paper
precmat.ml <- chol2inv(chol(dat.covmat))
#precmat.ml = (Var(Y|all parameters) + H%*%var.beta%*%t(H))^-1
rm(dat.covmat)

H.c <- as.matrix(m.h.c.rev)
H.e <- as.matrix(m.h.e.rev)
H.ml.lambda <- as.matrix(v.h.rev)
H <- as.matrix(m.h.test)

var.bc <- var.beta[1:ncol(m.h.c),1:ncol(m.h.c)]
var.be <- var.beta[(ncol(m.h.c)+1):(ncol(m.h)), (ncol(m.h.c)+1):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
covmat.ml.lambda <- H.ml.lambda %*% var.beta[1:ncol(m.h.e),1:ncol(m.h.e)] %*% t(H.ml.lambda) + cov.x1.x2(x.e.std.rev, x.e.std.rev, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.e.std.rev)[1])
precmat.ml.lambda <- chol2inv(chol(covmat.ml.lambda))
crosscov.ml.lambda <- gp.crosscov(x.v2, x.e.std.rev, pars$v_sigma, pars$v_theta, H, H.ml.lambda, var.beta[1:ncol(m.h.e),1:ncol(m.h.e)])
lambda.ml <- exp( crosscov.ml.lambda %*%( precmat.ml.lambda %*% pars$logLambda) )
var.ml.full <- var.ml + diag(lambda.ml)

MSE.mlCV3 <- MSE(y.test, mean.ml)

x.test.hetmlCV3<- x.test
y.test.hetmlCV3<- y.test

mean.mlCV3= mean.ml
lambda.mlCV3= lambda.ml

save(MSE.mlCV3, mean.mlCV3, lambda.mlCV3, x.test.hetmlCV3, y.test.hetmlCV3,
     file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV3hetml_calculated_quantities.RData")


################################################################################
################################################################################

#4th test set
test.inds<- randOrderInds[61:80]

################################################################################
#x.test<- rbind(x.true,parVals10m[test.inds,-1])
x.test<- rbind(parVals10m[test.inds,-1])
x.hetgp<- parVals10m[-test.inds,-1]
x.c <- parVals50m[-test.inds,-1]
x.e <- parVals10m[-test.inds,-1]

# load data
#y data
#y.test1 <- c(0,metrics.10m$EuclideanDists[test.inds]) #true Euclidean distance between output and reality
y.test1 <- metrics.10m$EuclideanDists[test.inds]
y.hetgp1<- metrics.10m$EuclideanDists[-test.inds]
y.c1<- metrics.10mfrom50m$EuclideanDists[-test.inds]
y.e1<- metrics.10m$EuclideanDists[-test.inds]
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
n.c <- length(y.c2)
n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

#create H
m.h.c.rev<- as.matrix(m.h.c[n.c:1,])

m.h.0<- matrix(0, ncol = ncol(m.h.c.rev), nrow = nrow(m.h.c.rev) - n.e)
m.h.e.rev<- as.matrix(m.h.e[n.e:1,])
m.h.2<- rbind(m.h.0, m.h.e.rev)

m.h <- cbind(m.h.c.rev, m.h.2)

v.h.rev<- cbind(1, x.e.std.rev)

#tail(m.h)
var.beta <- diag(1, ncol(m.h))
var.beta.m.h.c<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]
var.beta.m.h.v<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]

ml.data <- list(
  ## data ##
  m_p = ncol(m.h.e), 
  v_p = ncol(v.h.rev),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h, 
  v_H = v.h.rev,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  m_beta_m = rep(0, ncol(m.h.e)), 
  m_beta_s = var.beta.m.h.c,
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, 
  m_b_sigma = 2,
  m_nugget = 0,
  
  v_beta_m = rep(0, ncol(m.h.e)),
  v_beta_s = var.beta.m.h.v,
  v_a_theta = rep(2, ncol(x.test)), 
  v_b_theta = rep(1, ncol(x.test)),
  v_a_sigma = 2, 
  v_b_sigma = 2,
  v_nugget_a = 2, 
  v_nugget_b = 2,
  
  c_beta_m = rep(0, ncol(m.h.c)), 
  c_beta_s = var.beta.m.h.c,
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
  rstan::optimizing(sml, data = ml.data, verbose = F, as_vector = F)
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

save(pars,file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV4pars.hetML.RData")

load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV4pars.hetML.RData")

designmat <- matrix(0, ncol=ncol(m.h), nrow=n.c+n.e)
designmat[1:length(y.c),1:(ncol(m.h.c.rev))] <- m.h.c.rev
designmat[-(1:length(y.c)), (1:(ncol(m.h.e.rev)))] <- pars$rho*m.h.e.rev
designmat[-(1:length(y.c)), -(1:(ncol(m.h.e.rev)))] <- m.h.e.rev
#designmat is what is referred to as H on page 11 of the revised SML paper

dat.covmat <- dataCovmat(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
#dat.covmat = Var(Y|all parameters) + H%*%var.beta%*%t(H) on page 10 of revised SML paper
precmat.ml <- chol2inv(chol(dat.covmat))
#precmat.ml = (Var(Y|all parameters) + H%*%var.beta%*%t(H))^-1
rm(dat.covmat)

H.c <- as.matrix(m.h.c.rev)
H.e <- as.matrix(m.h.e.rev)
H.ml.lambda <- as.matrix(v.h.rev)
H <- as.matrix(m.h.test)

var.bc <- var.beta[1:ncol(m.h.c),1:ncol(m.h.c)]
var.be <- var.beta[(ncol(m.h.c)+1):(ncol(m.h)), (ncol(m.h.c)+1):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
covmat.ml.lambda <- H.ml.lambda %*% var.beta[1:ncol(m.h.e),1:ncol(m.h.e)] %*% t(H.ml.lambda) + cov.x1.x2(x.e.std.rev, x.e.std.rev, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.e.std.rev)[1])
precmat.ml.lambda <- chol2inv(chol(covmat.ml.lambda))
crosscov.ml.lambda <- gp.crosscov(x.v2, x.e.std.rev, pars$v_sigma, pars$v_theta, H, H.ml.lambda, var.beta[1:ncol(m.h.e),1:ncol(m.h.e)])
lambda.ml <- exp( crosscov.ml.lambda %*%( precmat.ml.lambda %*% pars$logLambda) )
var.ml.full <- var.ml + diag(lambda.ml)

MSE.mlCV4 <- MSE(y.test, mean.ml)

x.test.hetmlCV4<- x.test
y.test.hetmlCV4<- y.test

mean.mlCV4= mean.ml
lambda.mlCV4= lambda.ml

save(MSE.mlCV4, mean.mlCV4, lambda.mlCV4, x.test.hetmlCV4, y.test.hetmlCV4,
     file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV4hetml_calculated_quantities.RData")

################################################################################
################################################################################

#5th test set
test.inds<- randOrderInds[81:100]

################################################################################
#x.test<- rbind(x.true,parVals10m[test.inds,-1])
x.test<- rbind(parVals10m[test.inds,-1])
x.hetgp<- parVals10m[-test.inds,-1]
x.c <- parVals50m[-test.inds,-1]
x.e <- parVals10m[-test.inds,-1]

# load data
#y data
#y.test1 <- c(0,metrics.10m$EuclideanDists[test.inds]) #true Euclidean distance between output and reality
y.test1 <- metrics.10m$EuclideanDists[test.inds]
y.hetgp1<- metrics.10m$EuclideanDists[-test.inds]
y.c1<- metrics.10mfrom50m$EuclideanDists[-test.inds]
y.e1<- metrics.10m$EuclideanDists[-test.inds]
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
n.c <- length(y.c2)
n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

#create H
m.h.c.rev<- as.matrix(m.h.c[n.c:1,])

m.h.0<- matrix(0, ncol = ncol(m.h.c.rev), nrow = nrow(m.h.c.rev) - n.e)
m.h.e.rev<- as.matrix(m.h.e[n.e:1,])
m.h.2<- rbind(m.h.0, m.h.e.rev)

m.h <- cbind(m.h.c.rev, m.h.2)

v.h.rev<- cbind(1, x.e.std.rev)

#tail(m.h)
var.beta <- diag(1, ncol(m.h))
var.beta.m.h.c<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]
var.beta.m.h.v<- var.beta[1:ncol(m.h.e),1:ncol(m.h.e)]

ml.data <- list(
  ## data ##
  m_p = ncol(m.h.e), 
  v_p = ncol(v.h.rev),
  N = n.e + n.c, 
  K = ncol(x.test),
  n_c = n.c, 
  n_e = n.e,
  x_e = x.e.std.rev, 
  x_c = x.c.std.rev,
  m_H = m.h, 
  v_H = v.h.rev,
  y_c = y.c2, 
  y_e = y.e2,
  
  ## priors ##
  m_beta_m = rep(0, ncol(m.h.e)), 
  m_beta_s = var.beta.m.h.c,
  m_a_theta = rep(2, ncol(x.test)), 
  m_b_theta = rep(1, ncol(x.test)),
  m_a_sigma = 2, 
  m_b_sigma = 2,
  m_nugget = 0,
  
  v_beta_m = rep(0, ncol(m.h.e)),
  v_beta_s = var.beta.m.h.v,
  v_a_theta = rep(2, ncol(x.test)), 
  v_b_theta = rep(1, ncol(x.test)),
  v_a_sigma = 2, 
  v_b_sigma = 2,
  v_nugget_a = 2, 
  v_nugget_b = 2,
  
  c_beta_m = rep(0, ncol(m.h.c)), 
  c_beta_s = var.beta.m.h.c,
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
  rstan::optimizing(sml, data = ml.data, verbose = F, as_vector = F)
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

save(pars,file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV5pars.hetML.RData")

load(".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV5pars.hetML.RData")

designmat <- matrix(0, ncol=ncol(m.h), nrow=n.c+n.e)
designmat[1:length(y.c),1:(ncol(m.h.c.rev))] <- m.h.c.rev
designmat[-(1:length(y.c)), (1:(ncol(m.h.e.rev)))] <- pars$rho*m.h.e.rev
designmat[-(1:length(y.c)), -(1:(ncol(m.h.e.rev)))] <- m.h.e.rev
#designmat is what is referred to as H on page 11 of the revised SML paper

dat.covmat <- dataCovmat(x.c.std.rev, x.e.std.rev, pars, var.beta,  designmat)
#dat.covmat = Var(Y|all parameters) + H%*%var.beta%*%t(H) on page 10 of revised SML paper
precmat.ml <- chol2inv(chol(dat.covmat))
#precmat.ml = (Var(Y|all parameters) + H%*%var.beta%*%t(H))^-1
rm(dat.covmat)

H.c <- as.matrix(m.h.c.rev)
H.e <- as.matrix(m.h.e.rev)
H.ml.lambda <- as.matrix(v.h.rev)
H <- as.matrix(m.h.test)

var.bc <- var.beta[1:ncol(m.h.c),1:ncol(m.h.c)]
var.be <- var.beta[(ncol(m.h.c)+1):(ncol(m.h)), (ncol(m.h.c)+1):(ncol(m.h))]
crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be)
#assume a priori E(y|x) = 0 for any x
mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))

var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,2)
var.ml <- var.ml + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
covmat.ml.lambda <- H.ml.lambda %*% var.beta[1:ncol(m.h.e),1:ncol(m.h.e)] %*% t(H.ml.lambda) + cov.x1.x2(x.e.std.rev, x.e.std.rev, pars$v_sigma, pars$v_theta, 2) + diag(pars$v_nugget, dim(x.e.std.rev)[1])
precmat.ml.lambda <- chol2inv(chol(covmat.ml.lambda))
crosscov.ml.lambda <- gp.crosscov(x.v2, x.e.std.rev, pars$v_sigma, pars$v_theta, H, H.ml.lambda, var.beta[1:ncol(m.h.e),1:ncol(m.h.e)])
lambda.ml <- exp( crosscov.ml.lambda %*%( precmat.ml.lambda %*% pars$logLambda) )
var.ml.full <- var.ml + diag(lambda.ml)

MSE.mlCV5 <- MSE(y.test, mean.ml)

x.test.hetmlCV5<- x.test
y.test.hetmlCV5<- y.test

mean.mlCV5= mean.ml
lambda.mlCV5= lambda.ml

save(MSE.mlCV5, mean.mlCV5, lambda.mlCV5, x.test.hetmlCV5, y.test.hetmlCV5,
     file=".../FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV5hetml_calculated_quantities.RData")


################################################################################
################################################################################
