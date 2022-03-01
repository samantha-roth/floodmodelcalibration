library("foreach")
library("doParallel")

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha")
wd <- getwd()

#sample_river <-
#  as.matrix(read.delim2("./LISFLOOD/Sample_Selinsgrove2011.river", header = F)) # sample Susquehanna River file
sample_river <-
  as.matrix(read.delim2("./LISFLOOD/Sample_Selinsgrove.river", header = F)) # sample Susquehanna River file

sample_par <-
  as.matrix(read.delim2("./LISFLOOD/Sample_Selinsgrove10.par", header = F)) # sample parameter file



load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_thin500.RData")

# model setup and run
DEM_file <- "dem10.asc"

emulator<- "homML"

if(emulator=="hetML"){
  mat<- res_hetML_thin
}

if(emulator=="homML"){
  mat<- res_homML_thin
}

if(emulator=="hetGP"){
  mat<- res_hetGP_thin
}

if(emulator=="homGP"){
  mat<- res_homGP_thin
}

if(emulator=="hetGPCheap"){
  mat<- res_hetGPCheap_thin
}

if(emulator=="homGPCheap"){
  mat<- res_homGPCheap_thin
}



colnames(mat)[1:2]<-c("n_ch","n_fp")
#parameter sets and functions
source('./flood_extent_function_4pars_10m_MCMCpars.R')


# Output folders
if(emulator=="hetML"){
  
  if (dir.exists(paste0(wd, "/Outputs10m")) == F)
    dir.create(paste0(wd, "/Outputs10m"))
  if (dir.exists(paste0(wd, "/Outputs10m/hetML")) == F)
    dir.create(paste0(wd, "/Outputs10m/hetML"))
  if (dir.exists(paste0(wd, "/Outputs10m/hetML/Extent")) == F)
    dir.create(paste0(wd, "/Outputs10m/hetML/Extent"))
}

if(emulator=="homML"){
  
  if (dir.exists(paste0(wd, "/Outputs10m")) == F)
    dir.create(paste0(wd, "/Outputs10m"))
  if (dir.exists(paste0(wd, "/Outputs10m/homML")) == F)
    dir.create(paste0(wd, "/Outputs10m/homML"))
  if (dir.exists(paste0(wd, "/Outputs10m/homML/Extent")) == F)
    dir.create(paste0(wd, "/Outputs10m/homML/Extent"))
}

if(emulator=="hetGP"){
  
  if (dir.exists(paste0(wd, "/Outputs10m")) == F)
    dir.create(paste0(wd, "/Outputs10m"))
  if (dir.exists(paste0(wd, "/Outputs10m/hetGP")) == F)
    dir.create(paste0(wd, "/Outputs10m/hetGP"))
  if (dir.exists(paste0(wd, "/Outputs10m/hetGP/Extent")) == F)
    dir.create(paste0(wd, "/Outputs10m/hetGP/Extent"))
}

if(emulator=="homGP"){
  
  if (dir.exists(paste0(wd, "/Outputs10m")) == F)
    dir.create(paste0(wd, "/Outputs10m"))
  if (dir.exists(paste0(wd, "/Outputs10m/homGP")) == F)
    dir.create(paste0(wd, "/Outputs10m/homGP"))
  if (dir.exists(paste0(wd, "/Outputs10m/homGP/Extent")) == F)
    dir.create(paste0(wd, "/Outputs10m/homGP/Extent"))
}

if(emulator=="hetGPCheap"){
  
  if (dir.exists(paste0(wd, "/Outputs10m")) == F)
    dir.create(paste0(wd, "/Outputs10m"))
  if (dir.exists(paste0(wd, "/Outputs10m/hetGPCheap")) == F)
    dir.create(paste0(wd, "/Outputs10m/hetGPCheap"))
  if (dir.exists(paste0(wd, "/Outputs10m/hetGPCheap/Extent")) == F)
    dir.create(paste0(wd, "/Outputs10m/hetGPCheap/Extent"))
}

if(emulator=="homGPCheap"){
  
  if (dir.exists(paste0(wd, "/Outputs10m")) == F)
    dir.create(paste0(wd, "/Outputs10m"))
  if (dir.exists(paste0(wd, "/Outputs10m/homGPCheap")) == F)
    dir.create(paste0(wd, "/Outputs10m/homGPCheap"))
  if (dir.exists(paste0(wd, "/Outputs10m/homGPCheap/Extent")) == F)
    dir.create(paste0(wd, "/Outputs10m/homGPCheap/Extent"))
}



run_start = 1 #starting row number of the parameters table to read

run_end = 500 #ending row number of the parameters table to read



#setwd(paste0(wd,'/LISFLOOD'))
#system("./chmod a+x lisflood.exe") # activate the .exe file for the first run

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) # -1 not to overload system
registerDoParallel(cl)

start<-proc.time()
foreach (i = run_start:run_end) %do% {
  print(paste("Run", i))
  flood_extent_run(i)
}



stopCluster(cl)
end<-proc.time()

print(end-start)
