# Set working directory
#(wd <- getwd())
#if (!is.null(wd))
#  setwd(wd)


library(raster)

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha")
wd <- getwd()

sample_river <-
  as.matrix(read.delim2("./LISFLOOD/Sample_Selinsgrove2011truth.river", header = F)) # sample Susquehanna River file
sample_par <-
  as.matrix(read.delim2("./LISFLOOD/Sample_Selinsgrove50truth.par", header = F)) # sample parameter file

# model setup and run
DEM_file <- "dem50.asc" # 1m digital elevation model
mat<-data.frame(0.0305,0.045) #assuming a single value of 0.03 for both river bed and floodplain roughness coefficients
colnames(mat)<-c("n_ch","n_fp")

# calling the flood extent function
source('./flood_extent_function_50mtruth.R')

# Output folders
if (dir.exists(paste0(wd, "/Outputs50m")) == F)
  dir.create(paste0(wd, "/Outputs50m"))
if (dir.exists(paste0(wd, "/Outputs50m/Extent")) == F)
  dir.create(paste0(wd, "/Outputs50m/Extent"))

i=1 #this is a single run with the 50m dem

#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/LISFLOOD")
#system("./chmod a+x lisflood.exe") # activate the .exe file for the first run

#setwd(paste0(wd,'/LISFLOOD'))
#system("./chmod a+x lisflood.exe") # activate the .exe file for the first run
#system("chmod a+x lisflood.exe") # activate the .exe file for the first run


#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha")
#wd<- getwd()

flood_extent_run(i)
