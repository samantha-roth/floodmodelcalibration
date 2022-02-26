#Compare single resolution emulation results 

rm(list=ls())
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV1homGP_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2homGP_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV3homGP_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV4homGP_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV5homGP_calculated_quantities.RData")

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV1hetGP_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2hetGP_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV3hetGP_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV4hetGP_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV5hetGP_calculated_quantities.RData")


mean(c(sqrt(MSE.hetgpCV1),sqrt(MSE.hetgpCV2),sqrt(MSE.hetgpCV3),sqrt(MSE.hetgpCV4),sqrt(MSE.hetgpCV5)))
mean(c(sqrt(MSE.gpCV1),sqrt(MSE.gpCV2),sqrt(MSE.gpCV3),sqrt(MSE.gpCV4),sqrt(MSE.gpCV5)))

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV1hetml_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV2hetml_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV3hetml_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV4hetml_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/hold20/prior1/10mfrom50m/CV5hetml_calculated_quantities.RData")

mean(c(sqrt(MSE.mlCV1),sqrt(MSE.mlCV2),sqrt(MSE.mlCV3),sqrt(MSE.mlCV4),sqrt(MSE.mlCV5)))

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV1homML_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV2homML_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV3homML_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV4homML_calculated_quantities.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/EuclideanDistance/holdTruePlus10/prior1/10mfrom50m/CV5homML_calculated_quantities.RData")

mean(c(sqrt(MSE.hommlCV1),sqrt(MSE.hommlCV2),sqrt(MSE.hommlCV3),sqrt(MSE.hommlCV4),sqrt(MSE.hommlCV5)))

