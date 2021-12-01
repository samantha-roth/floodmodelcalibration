#My Heteroscedastic Gaussian Process Emulator
rm(list=ls())
library(raster)
#load the data

load("C:/FloodingModelCalibrationProject/ensemble_uniform.RData")
Manning<- mat; rm(mat) #different values for each of the 100 runs of Manning's
#roughness coefficient
load("C:/FloodingModelCalibrationProject/data/Extent/allRunArray.RData")

# the dimension of allRunArray[,,i] is number x grid cells by number y grid cells

Run_1= raster("C:/FloodingModelCalibrationProject/data/Extent/Run_1.asc")
plot(Run_1)
Run1Mat<- as.matrix(Run_1)

#get the coordinates of cells
coords <- xyFromCell(Run_1,1:ncell(Run_1))
x<- rev(unique(coords[,2])) #longitude
y<- unique(coords[,1]) #latitude

vals1<- extract(Run_1,coords)
yIndsIwant1<- which(coords[,1]>343100)
yIndsIwant2<- which(coords[,1]<343600)
xIndsIwant1<- which(coords[,2]>4519000)
xIndsIwant2<- which(coords[,2]<4519200)
yIndsIwant<- intersect(yIndsIwant1,yIndsIwant2)
xIndsIwant<- intersect(xIndsIwant1,xIndsIwant2)
coordsIwantInds<- intersect(yIndsIwant,xIndsIwant)
coordsIwant<- coords[coordsIwantInds,]
x.s<- rev(unique(coordsIwant[,2])) #longitude
y.s<- unique(coordsIwant[,1]) #latitude
Run_1smaller<- rasterFromCells(Run_1,coordsIwantInds)
values(Run_1smaller)<- vals1[coordsIwantInds]
plot(Run_1smaller)

#r <- raster(ncols=100, nrows=100)
#values(r)<- rnorm(100^2)
#plot(r)

#coordsEx<- xyFromCell(r,1:ncell(r))

#cells <- which(coordsEx[,1]<50)
#r2 <- rasterFromCells(r, cells, values=FALSE)
#valsr<- extract(r,xyFromCell(r,1:ncell(r)))
#copyvals<- valsr[cells]
#values(r2)<- copyvals
#plot(r2)



#valsr2<- extract(r2,xyFromCell(r2,1:ncell(r2)))
#cbind(1:ncell(r2), getValues(r2))


################################################################################
#DATA FORMATTING FOR USING IN HETGP MODEL#

nRuns<- 100
nx<- length(x)
ny<- length(y)

allY<- rep(coords[,1], nRuns)
allX<- rep(coords[,2], nRuns)
allvals<- rep(NA, nx*ny*nRuns)
runInd<- rep(1:nRuns,each=nx*ny)
manningInd<- rep(Manning,each=nx*ny)

for(i in 1:100){
  run<- raster(paste("C:/FloodingModelCalibrationProject/data/Extent/Run_",i,".asc",sep=""))
  vals<- extract(run,coords)
  allvals[((i-1)*nx*ny+1):(i*nx*ny)]<- vals
}

allRunData<- data.frame(
  "run"= runInd,
  "manning"= manningInd,
  "val"= allvals,
  "y"= allY,
  "x"= allX
)

save(allRunData,file="C:/FloodingModelCalibrationProject/data/allRunData.RData")
################################################################################
#smaller version of data above

nRuns<- 100
nx.s<- length(x.s)
ny.s<- length(y.s)

allY<- rep(coordsIwant[,1], nRuns)
allX<- rep(coordsIwant[,2], nRuns)
allvals<- rep(NA, nx.s*ny.s*nRuns)
runInd<- rep(1:nRuns,each=nx.s*ny.s)
manningInd<- rep(Manning,each=nx.s*ny.s)

for(i in 1:100){
  run<- raster(paste("C:/FloodingModelCalibrationProject/data/Extent/Run_",i,".asc",sep=""))
  vals1<- extract(run,coords)
  run_smaller<- rasterFromCells(run,coordsIwantInds)
  values(run_smaller)<- vals1[coordsIwantInds]
  vals<- extract(run_smaller,coordsIwant)
  allvals[((i-1)*nx.s*ny.s+1):(i*nx.s*ny.s)]<- vals
}

allRunDataS<- data.frame(
  "run"= runInd,
  "manning"= manningInd,
  "val"= allvals,
  "y"= allY,
  "x"= allX
)

save(allRunDataS,nx.s,ny.s,nRuns,file="C:/FloodingModelCalibrationProject/data/allRunDataS.RData")
