#fit HetGP and HetML to 4 parameter sample
rm(list=ls())
setwd("C:/FloodingModelCalibrationProject/sml-athena-main")
library(rstan)
library(boot)
library("R.matlab")
source("GPfunctionsOptim.R")
source("hetGPfunctions.R")
hetGP <- stan_model("hetGP.stan")
sml <- stan_model("het-SML.stan")

#load parameters
load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/allParVals.RData")
parVals10m<- as.data.frame(parVals)
load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/allParVals.RData")
parVals50m<- as.data.frame(parVals)


#load predictions
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10m.RData")
#50m disaggregated
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10mfrom50m.RData")
#10m aggregated
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics50mfrom10m.RData")


#load true values
#true parameter values
parsTrue<-read.csv("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/parVals/RunTrue_1.csv")


nPars=4

#x data
x.true<- parsTrue[,-(1:2)]

#randomly hold out data to test

################################################################################
#for when testing data is used
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

################################################################################
#no heldout data

#x.test<- as.matrix(x.true)
#x.hetgp<- parVals10m[,-1]
#x.c <- parVals50m[,-1]
#x.e <- parVals10m[,-1]

# load data
#y data
#y.test1 <- 0 #true Euclidean distance between output and reality
#y.hetgp1<- metrics.10m$EuclideanDists
#50m disaggregated
#y.c1<- metrics.10mfrom50m$EuclideanDists
#10m aggregated
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

#save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/pars.hetGP.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/pars.hetGP.RData")

#save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/pars.hetGP.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/pars.hetGP.RData")

#save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/pars.hetGP.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/pars.hetGP.RData")


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
plot(lambda.hetgp, y.test)
mean((mean.hetgp - y.test)^2)
sqrt(mean((mean.hetgp - y.test)^2))

plot(x.test[,1],mean.hetgp); plot(x.test[,1],lambda.hetgp)

#mean of predictive distribution of emulator doesn't capture the increase in ED at lowest channel roughness values
#mean of predictive distribution of lambda^2_e is very large for lowest channel roughness values

MSE.hetgp <- MSE(y.test, mean.hetgp)


y.test.hetgp<- y.test
x.test.hetgp<- x.test


#save(mean.hetgp, lambda.hetgp, x.test.hetgp, y.test.hetgp, MSE.hetgp, file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/hetGP_calculated_quantities.RData")
#save(mean.hetgp, lambda.hetgp, x.test.hetgp, y.test.hetgp, MSE.hetgp, file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/hetGP_calculated_quantities.RData")

save(mean.hetgp, lambda.hetgp, x.test.hetgp, y.test.hetgp, MSE.hetgp,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/hetGP_calculated_quantities.RData")


load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/hetGP_calculated_quantities.RData")


################################################################################
########################### fit SML emulator ###################################
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
  
  
  #disaggregate
  #save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/pars.hetML.RData")
  #load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/pars.hetML.RData")
  
  #aggregate
  #save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/pars.hetML.RData")
  #load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/pars.hetML.RData")

  
save(pars,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/pars.hetML.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/pars.hetML.RData")


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

MSE.ml <- MSE(y.test, mean.ml)

x.test.hetml<- x.test
y.test.hetml<- y.test

save(MSE.ml,mean.ml,lambda.ml,x.test.hetml,y.test.hetml,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/hetml_calculated_quantities.RData")

#save(MSE.ml,mean.ml,lambda.ml,x.test.hetml,y.test.hetml,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/hetml_calculated_quantities.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/hetml_calculated_quantities.RData")


#save(MSE.ml,mean.ml, lambda.ml, x.test.hetml, y.test.hetml,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/hetml_calculated_quantities.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTrue/prior1/10mfrom50m/ag/hetml_calculated_quantities.RData")

#### compare the fits ####

## score & MSE
#mse
MSE.hetgp <- MSE(y.test, mean.hetgp)
MSE.sml <- MSE(y.test, mean.ml)
MSE.hetgp; MSE.sml
sqrt(MSE.hetgp); sqrt(MSE.sml)

#MAE
MAE.hetgp <- MAE(y.test, mean.hetgp)
MAE.sml <- MAE(y.test, mean.ml)
MAE.hetgp; MAE.sml

## Rmse
sqrt(MSE.hetgp); sqrt(MSE.sml)

#score
Score.hetgp <- Score(y.test, mean.hetgp, diag(var.hetgp.full))
Score.sml <- Score(y.test, mean.ml, diag(var.ml.full))
sum(Score.hetgp); sum(Score.sml)
sumScore.hetgp<- sum(Score.hetgp)
sumScore.sml<- sum(Score.sml)

#save values of metrics
save(MSE.hetgp,MSE.sml,mean.hetgp,mean.ml,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/het_performance.RData")

#plot fits
plot(x.test[,1], mean.ml, main="E(Y_test|Y) from HetSML Emulator",xlab="Channel roughness",ylab="Predicted Euclidean Distance")
plot(x.test[,2], mean.ml, main="E(Y_test|Y) from HetSML Emulator",xlab="Floodplain roughness",ylab="Predicted Euclidean Distance")
plot(x.test[,3], mean.ml, main="E(Y_test|Y) from HetSML Emulator",xlab="River width error",ylab="Predicted Euclidean Distance")
plot(x.test[,4], mean.ml, main="E(Y_test|Y) from HetSML Emulator",xlab="Riverbed elevation error",ylab="Predicted Euclidean Distance")
plot(mean.ml, y.test, main="E(Y_test|Y) from HetSML Emulator",xlab="Predicted Euclidean Distance",ylab="True Euclidean Distance")

plot(x.test[,1], mean.hetgp, main="E(Y_test|Y) from HetGP Emulator",xlab="Channel roughness",ylab="Predicted Euclidean Distance")
plot(x.test[,2], mean.hetgp, main="E(Y_test|Y) from HetGP Emulator",xlab="Floodplain roughness",ylab="Predicted Euclidean Distance")
plot(x.test[,3], mean.hetgp, main="E(Y_test|Y) from HetGP Emulator",xlab="River width error",ylab="Predicted Euclidean Distance")
plot(x.test[,4], mean.hetgp, main="E(Y_test|Y) from HetGP Emulator",xlab="Riverbed elevation error",ylab="Predicted Euclidean Distance")
plot(mean.hetgp, y.test, main="E(Y_test|Y) from HetGP Emulator",xlab="Predicted Euclidean Distance",ylab="True Euclidean Distance")

#plot predicted lambda^2_e values

## look at residuals
par(mfrow=c(1,2))

chol.het <- forwardsolve(t(chol(var.hetgp.full)), (y.test - mean.hetgp))
chol.ml <- forwardsolve(t(chol(var.ml.full)), (y.test - mean.ml))

## coverage plot


coverage <- function(resids){
  alpha <- seq(0, 1, length = length(resids))
  emp.cov <- rep(0, length(resids))
  for(i in 1:length(resids)){
    emp.cov[i] <- sum(abs(resids) < qnorm(1-alpha[i]/2))/length(resids)
  }
  list(y=emp.cov,x = 1 - alpha) 
}


