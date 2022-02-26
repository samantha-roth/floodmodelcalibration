# floodmodelcalibration

To calculate Euclidean Distances between model output and truth:
- Run metricEvaluation10m50m_4pars.R

To emulate:
- Fit HomMR and HomGP10 emulators: fitHomML.R 
- Fit HetMR and HetGP10 emulators: fitHetML.R
- Fit HomGP50 emulator: fitHomGPCheap.R
- Fit HetGP50 emulator fitHetGPCheap.R

Emulator cross-validation study:
- HomMR: fitHomML_CVallpars.R
- HetMR: fitHetML_CVallpars.R
- HetGP10: fitHetGP_CV1.R, fitHetGP_CV2.R, fitHetGP_CV3.R, fitHetGP_CV4.R, fitHetGP_CV5.R
- HomGP10: fitHomGP_CV1.R, fitHomGP_CV2.R, fitHomGP_CV3.R, fitHomGP_CV4.R, fitHomGP_CV5.R
- HomGP50: fitHomGPCheap_CVallParallel.R (function to source: homGPCheapEmulatorCVPredsFunction.R)
- HetGP50: fitHetGPCheap_CVallParallel.R (function to source: hetGPCheapEmulatorCVPredsFunction.R)
- Compare emulators for 10m model: compareEmulationHoldOut20.R


To calibrate:
- Calibrate using HomMR: homML_cal.R
- Calibrate using HomGP10: homGP_cal.R
- Calibrate using HomGP50: homGPCheap_cal.R

Compare calibration results:
- Look at ESS, ACF plots, density plots of first and second half of Markov Chains, posterior densities: 
  - compareCalibrationResults.R 
  - plotPosteriorDensities.R
- Look at joint posterior densities: plotJointPosteriorDensities.R

Predictions for 2011:


Predictions for 2004:
