library("foreach")
library("doParallel")

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha")
wd <- getwd()

sample_river <-
  as.matrix(read.delim2("./LISFLOOD/Sample_Selinsgrove2011.river", header = F)) # sample Susquehanna River file
sample_par <-
  as.matrix(read.delim2("./LISFLOOD/Sample_Selinsgrove10.par", header = F)) # sample parameter file

prior= 1
set= "second"

#load('./lhs_samples.RData')
#load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/4param_lhs_samples.RData")
#load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/4param_lhs_samples2.RData")
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/4param_lhs_samples_allU.RData")
#load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/4param_lhs_samples_tN.RData")
#load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/4param_lhs_samples_allU_Extra.RData")


# model setup and run
DEM_file <- "dem10.asc"
if(set=="first") mat<-samp.E
if(set=="second") mat<-samp.C[201:300,]
if(set=="third") mat<-samp.C[301:400,]
if(set=="fourth") mat<-samp.C[401:500,]
if(set=="fifth") mat<-samp.C[501:600,]
if(set=="sixth") mat<-samp.C[601:700,]
if(set=="seventh") mat<-samp.C[701:800,]
colnames(mat)[1:2]<-c("n_ch","n_fp")

#parameter sets and functions
source('./flood_extent_function_4pars_10m.R')


# Output folders
if(prior==1){

if (dir.exists(paste0(wd, "/Outputs10m")) == F)
  dir.create(paste0(wd, "/Outputs10m"))
if (dir.exists(paste0(wd, "/Outputs10m/prior1")) == F)
  dir.create(paste0(wd, "/Outputs10m/prior1"))
if (dir.exists(paste0(wd, "/Outputs10m/prior1/Extent")) == F)
  dir.create(paste0(wd, "/Outputs10m/prior1/Extent"))

	if(set=="second"){
  	if (dir.exists(paste0(wd, "/Outputs10m/prior1/second")) == F)
    		dir.create(paste0(wd, "/Outputs10m/prior1/second"))
  	if (dir.exists(paste0(wd, "/Outputs10m/prior1/second/Extent")) == F)
    		dir.create(paste0(wd, "/Outputs10m/prior1/second/Extent"))
  	}

	if(set=="third"){
  	if (dir.exists(paste0(wd, "/Outputs10m/prior1/third")) == F)
    		dir.create(paste0(wd, "/Outputs10m/prior1/third"))
  	if (dir.exists(paste0(wd, "/Outputs10m/prior1/third/Extent")) == F)
    		dir.create(paste0(wd, "/Outputs10m/prior1/third/Extent"))
  	}
	
	if(set=="fourth"){
  	if (dir.exists(paste0(wd, "/Outputs10m/prior1/fourth")) == F)
    		dir.create(paste0(wd, "/Outputs10m/prior1/fourth"))
  	if (dir.exists(paste0(wd, "/Outputs10m/prior1/fourth/Extent")) == F)
    		dir.create(paste0(wd, "/Outputs10m/prior1/fourth/Extent"))
	}
  
  if(set=="fifth"){
    if (dir.exists(paste0(wd, "/Outputs10m/prior1/fifth")) == F)
      dir.create(paste0(wd, "/Outputs10m/prior1/fifth"))
    if (dir.exists(paste0(wd, "/Outputs10m/prior1/fifth/Extent")) == F)
      dir.create(paste0(wd, "/Outputs10m/prior1/fifth/Extent"))
  }
  
  if(set=="sixth"){
    if (dir.exists(paste0(wd, "/Outputs10m/prior1/sixth")) == F)
      dir.create(paste0(wd, "/Outputs10m/prior1/sixth"))
    if (dir.exists(paste0(wd, "/Outputs10m/prior1/sixth/Extent")) == F)
      dir.create(paste0(wd, "/Outputs10m/prior1/sixth/Extent"))
  }
  
  if(set=="seventh"){
    if (dir.exists(paste0(wd, "/Outputs10m/prior1/seventh")) == F)
      dir.create(paste0(wd, "/Outputs10m/prior1/seventh"))
    if (dir.exists(paste0(wd, "/Outputs10m/prior1/seventh/Extent")) == F)
      dir.create(paste0(wd, "/Outputs10m/prior1/seventh/Extent"))
  }
  
}

if(prior==2){

if (dir.exists(paste0(wd, "/Outputs10m")) == F)
  dir.create(paste0(wd, "/Outputs10m"))
if (dir.exists(paste0(wd, "/Outputs10m/prior2")) == F)
  dir.create(paste0(wd, "/Outputs10m/prior2"))
if (dir.exists(paste0(wd, "/Outputs10m/prior2/Extent")) == F)
  dir.create(paste0(wd, "/Outputs10m/prior2/Extent"))
}



run_start = 1 #starting row number of the parameters table to read
run_end = nrow(mat) #ending row number of the parameters table to read

#run_start = nrow(mat)/2 +1 #starting row number of the parameters table to read
#run_end = nrow(mat) #ending row number of the parameters table to read



#setwd(paste0(wd,'/LISFLOOD'))
#system("./chmod a+x lisflood.exe") # activate the .exe file for the first run

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) # -1 not to overload system
registerDoParallel(cl)

start<-proc.time()
#foreach (i = run_start:run_end) %do% {
#   print(paste("Run", i))
#   flood_extent_run(i)
#}

foreach (i = 6:7) %do% {
   print(paste("Run", i))
   flood_extent_run(i)
}


stopCluster(cl)
end<-proc.time()

print(end-start)

