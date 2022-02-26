#plot Lisflood calibrated preds

rm(list=ls())
#plot Lisflood calibration preds
library(raster)
library(RColorBrewer)

RunTrue.10m= raster("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/Extent/RunTrue_1.asc")
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
truevals.10m<- raster::extract(RunTrue.10m,coords.10m)
########################################################################################################################

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homGP.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetGP.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homML.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetML.RData")

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homGPCheap.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetGPCheap.RData")

########################################################################################################################


#Run_hetML_mean<- RunTrue.10m
#values(Run_hetML_mean)<- hetML_mean

Run_homML_mean<- RunTrue.10m
values(Run_homML_mean)<- homML_mean

#Run_hetGP_mean<- RunTrue.10m
#values(Run_hetGP_mean)<- hetGP_mean

Run_homGP_mean<- RunTrue.10m
values(Run_homGP_mean)<- homGP_mean


#Run_hetGPCheap_mean<- RunTrue.10m
#values(Run_hetGPCheap_mean)<- hetGPCheap_mean

Run_homGPCheap_mean<- RunTrue.10m
values(Run_homGPCheap_mean)<- homGPCheap_mean


#cuts=seq(0,9,by=.01)
#pal <- colorRampPalette(c("white","blue"))
#plot(RunTrue.10m, breaks=cuts, col = pal(length(cuts)+1)) #plot with defined breaks


#par(mfrow=c(3,2))
#plot(RunTrue.10m, main= "Observation", breaks=cuts, col = pal(length(cuts)+1))
#plot(Run_hetML_mean, main= "Mean prediction from HetMR emulation-calibration", breaks=cuts, col = pal(length(cuts)+1))
#plot(Run_homML_mean, main= "Mean prediction from HomMR emulation-calibration", breaks=cuts, col = pal(length(cuts)+1))
#plot(Run_homGP_mean, main= "Mean prediction from HomGP10 emulation-calibration", breaks=cuts, col = pal(length(cuts)+1))
#plot(Run_homGPCheap_mean, main= "Mean prediction from HomGP50 emulation-calibration", breaks=cuts, col = pal(length(cuts)+1))
#plot(Run_hetGPCheap_mean, main= "Mean prediction from HetGP50 emulation-calibration", breaks=cuts, col = pal(length(cuts)+1))



#plot the mean calibrated model runs with a map in the background
################################################################################

#GET MAP OF SELINSGROVE

library(ggmap)
library(terra)
library(raster)
library(rgdal)


RunTrue.10m= raster("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/Extent/RunTrue_1.asc")
#set crs
crs(RunTrue.10m)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
truevals.10m <- raster::extract(RunTrue.10m,coords.10m)

max.y.10m<- max(coords.10m[,2])
min.y.10m<- min(coords.10m[,2])
max.x.10m<- max(coords.10m[,1])
min.x.10m<- min(coords.10m[,1])

#Produce map of Selinsgrove
points<- coords.10m
v <- vect(points, crs="+proj=utm +zone=18 +datum=WGS84  +units=m")
v
y <- terra::project(v, "+proj=longlat +datum=WGS84")
y

lonlat <- geom(y)[, c("x", "y")]
head(lonlat, 3)

max.lat.10m<- max(lonlat[,2])
min.lat.10m<- min(lonlat[,2])
max.lon.10m<- max(lonlat[,1])
min.lon.10m<- min(lonlat[,1])


selingsgrove_box <- c(
  left = min.lon.10m,
  bottom = min.lat.10m,
  right = max.lon.10m,
  top = max.lat.10m
)

selinsgrove_stamen <- get_stamenmap(
  bbox = selingsgrove_box,
  maptype = "toner",
  zoom = 15
)

selinsgrove_ggmap<- ggmap(selinsgrove_stamen)

################################################################################
################################################################################
#good method,  now I just need to restrcit raster layer
#r <- raster(system.file("external/test.grd", package="raster"))


#plot true vals >0 

# just to make it reproducible with ggmap we have to transform to wgs84
RunTrue.10mG0<- RunTrue.10m
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
vals.10m<- raster::extract(RunTrue.10m,coords.10m)
ZeroInds<- which(vals.10m==0)
vals.10mG0<- vals.10m
vals.10mG0[ZeroInds]<- NA
values(RunTrue.10mG0)<- vals.10mG0
vals.10mincm<- as.integer(round(100*vals.10mG0))
values(RunTrue.10mG0)<- vals.10mincm #convert to cm since rtp can only deal with integer values

r <- projectRaster(RunTrue.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

dataType(r)
dataType(r)="INT4S"

rtp <- rasterToPolygons(r)

#selinsgrove_ggmap + 
#  geom_polygon(data = rtp, 
#               aes(x = long, y = lat, group = group, 
#                   fill = rep(rtp$RunTrue_1, each= 5)), 
#               #fill=layer),
#               size = 0, 
#               alpha = 0.5)  + 
#  scale_fill_gradientn("Flood Height (cm)", colors = heat.colors(500)) 


######This method works!!
################################################################################

#now do this for the mean predictions

################################################################################
#HetML
#crs(Run_hetML_mean)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
#Run_hetML_mean.10mG0<- Run_hetML_mean
#Run_hetML_meanvals.10m<- raster::extract(Run_hetML_mean,coords.10m)
#ZeroInds<- which(Run_hetML_meanvals.10m==0)
#Run_hetML_meanvals.10mG0<- Run_hetML_meanvals.10m
#Run_hetML_meanvals.10mG0[ZeroInds]<- NA
#values(Run_hetML_mean.10mG0)<- Run_hetML_meanvals.10mG0
#Run_hetML_meanvals.10mincm<- as.integer(round(100*Run_hetML_meanvals.10mG0))
#values(Run_hetML_mean.10mG0)<- Run_hetML_meanvals.10mincm #convert to cm since rtp can only deal with integer values

#r.hetML <- projectRaster(Run_hetML_mean.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
#dataType(r.hetML)="INT4S"

#rtp.hetML <- rasterToPolygons(r.hetML)

################################################################################
#HomML 

crs(Run_homML_mean)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homML_mean.10mG0<- Run_homML_mean
Run_homML_meanvals.10m<- raster::extract(Run_homML_mean,coords.10m)
ZeroInds<- which(Run_homML_meanvals.10m==0)
Run_homML_meanvals.10mG0<- Run_homML_meanvals.10m
Run_homML_meanvals.10mG0[ZeroInds]<- NA
values(Run_homML_mean.10mG0)<- Run_homML_meanvals.10mG0
Run_homML_meanvals.10mincm<- as.integer(round(100*Run_homML_meanvals.10mG0))
values(Run_homML_mean.10mG0)<- Run_homML_meanvals.10mincm #convert to cm since rtp can only deal with integer values

r.homML <- projectRaster(Run_homML_mean.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homML)="INT4S"

rtp.homML <- rasterToPolygons(r.homML)

################################################################################
#HetGP50

#crs(Run_hetGPCheap_mean)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
#Run_hetGPCheap_mean.10mG0<- Run_hetGPCheap_mean
#Run_hetGPCheap_meanvals.10m<- raster::extract(Run_hetGPCheap_mean,coords.10m)
#ZeroInds<- which(Run_hetGPCheap_meanvals.10m==0)
#Run_hetGPCheap_meanvals.10mG0<- Run_hetGPCheap_meanvals.10m
#Run_hetGPCheap_meanvals.10mG0[ZeroInds]<- NA
#values(Run_hetGPCheap_mean.10mG0)<- Run_hetGPCheap_meanvals.10mG0
#Run_hetGPCheap_meanvals.10mincm<- as.integer(round(100*Run_hetGPCheap_meanvals.10mG0))
#values(Run_hetGPCheap_mean.10mG0)<- Run_hetGPCheap_meanvals.10mincm #convert to cm since rtp can only deal with integer values

#r.hetGPCheap <- projectRaster(Run_hetGPCheap_mean.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
#dataType(r.hetGPCheap)="INT4S"

#rtp.hetGPCheap <- rasterToPolygons(r.hetGPCheap)


################################################################################
#HomGP50

crs(Run_homGPCheap_mean)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homGPCheap_mean.10mG0<- Run_homGPCheap_mean
Run_homGPCheap_meanvals.10m<- raster::extract(Run_homGPCheap_mean,coords.10m)
ZeroInds<- which(Run_homGPCheap_meanvals.10m==0)
Run_homGPCheap_meanvals.10mG0<- Run_homGPCheap_meanvals.10m
Run_homGPCheap_meanvals.10mG0[ZeroInds]<- NA
values(Run_homGPCheap_mean.10mG0)<- Run_homGPCheap_meanvals.10mG0
Run_homGPCheap_meanvals.10mincm<- as.integer(round(100*Run_homGPCheap_meanvals.10mG0))
values(Run_homGPCheap_mean.10mG0)<- Run_homGPCheap_meanvals.10mincm #convert to cm since rtp can only deal with integer values

r.homGPCheap <- projectRaster(Run_homGPCheap_mean.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homGPCheap)="INT4S"

rtp.homGPCheap <- rasterToPolygons(r.homGPCheap)


################################################################################
#HetGP10

#crs(Run_hetGP_mean)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
#Run_hetGP_mean.10mG0<- Run_hetGP_mean
#Run_hetGP_meanvals.10m<- raster::extract(Run_hetGP_mean,coords.10m)
#ZeroInds<- which(Run_hetGP_meanvals.10m==0)
#Run_hetGP_meanvals.10mG0<- Run_hetGP_meanvals.10m
#Run_hetGP_meanvals.10mG0[ZeroInds]<- NA
#values(Run_hetGP_mean.10mG0)<- Run_hetGP_meanvals.10mG0
#Run_hetGP_meanvals.10mincm<- as.integer(round(100*Run_hetGP_meanvals.10mG0))
#values(Run_hetGP_mean.10mG0)<- Run_hetGP_meanvals.10mincm #convert to cm since rtp can only deal with integer values

#r.hetGP <- projectRaster(Run_hetGP_mean.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
#dataType(r.hetGP)="INT4S"

#rtp.hetGP <- rasterToPolygons(r.hetGP)

################################################################################
#HomGP10

crs(Run_homGP_mean)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homGP_mean.10mG0<- Run_homGP_mean
Run_homGP_meanvals.10m<- raster::extract(Run_homGP_mean,coords.10m)
ZeroInds<- which(Run_homGP_meanvals.10m==0)
Run_homGP_meanvals.10mG0<- Run_homGP_meanvals.10m
Run_homGP_meanvals.10mG0[ZeroInds]<- NA
values(Run_homGP_mean.10mG0)<- Run_homGP_meanvals.10mG0
Run_homGP_meanvals.10mincm<- as.integer(round(100*Run_homGP_meanvals.10mG0))
values(Run_homGP_mean.10mG0)<- Run_homGP_meanvals.10mincm #convert to cm since rtp can only deal with integer values

r.homGP <- projectRaster(Run_homGP_mean.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homGP)="INT4S"

rtp.homGP <- rasterToPolygons(r.homGP)

################################################################################

#max(rtp.hetML$RunTrue_1)
max(rtp.homML$RunTrue_1)
#max(rtp.homGP$RunTrue_1)
max(rtp.homGPCheap$RunTrue_1)
#max(rtp.hetGPCheap$RunTrue_1)
max(rtp.homGP$RunTrue_1)

cuts=seq(0,9,by=.01)
pal <- colorRampPalette(c("white","blue"))
plot(RunTrue.10m, breaks=cuts, col = pal(length(cuts)+1)) #plot with defined breaks

#par(mfrow=c(3,2))
rtp$RunTrue_1<- rtp$RunTrue_1/100
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/Observation.jpeg",width = 800, height = 600)
selinsgrove_ggmap + 
  geom_polygon(data = rtp, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Flood Height (m)", colors = pal(900), limits= c(0,9)) +
  labs(title= "Observation") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()
beepr::beep()
################################################################################
rtp.homML$RunTrue_1<- rtp.homML$RunTrue_1/100
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomML_meanpred.jpeg",width = 800, height = 500)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homML, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homML$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Flood Height (m)", colors = pal(900), limits= c(0,9)) +
  labs(title= "Mean prediction from HomMR emulation-calibration")
dev.off()
beepr::beep()

################################################################################
rtp.homGPCheap$RunTrue_1<- rtp.homGPCheap$RunTrue_1/100
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomGP50_meanpred.jpeg",width = 800, height = 500)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homGPCheap, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homGPCheap$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Flood Height (m)", colors = pal(900), limits= c(0,9)) +
  labs(title= "Mean prediction from HomGP50 emulation-calibration") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24))
dev.off()
beepr::beep()
################################################################################
rtp.homGP$RunTrue_1<- rtp.homGP$RunTrue_1/100
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomGP10_meanpred.jpeg",width = 800, height = 500)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homGP, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homGP$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Flood Height (m)", colors = pal(900), limits= c(0,9)) +
  labs(title= "Mean prediction from HomGP10 emulation-calibration")+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24))
dev.off()
beepr::beep()
################################################################################
################################################################################
rtp.hetML$RunTrue_1<- rtp.hetML$RunTrue_1/100
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HetML_meanpred.jpeg",width = 800, height = 500)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.hetML, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.hetML$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900))+
  scale_fill_gradientn("Flood Height (m)", colors = pal(900), limits= c(0,9)) +
  labs(title= "Mean prediction from HetMR emulation-calibration")+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24))
dev.off()
beepr::beep()
################################################################################
rtp.hetGPCheap$RunTrue_1<- rtp.hetGPCheap$RunTrue_1/100
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HetGP50_meanpred.jpeg",width = 800, height = 500)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.hetGPCheap, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.hetGPCheap$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Flood Height (m)", colors = pal(900), limits= c(0,9)) +
  labs(title= "Mean prediction from HetGP50 emulation-calibration")+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24))
dev.off()
beepr::beep()
################################################################################

min(rtp.hetML$RunTrue_1)
min(rtp.homML$RunTrue_1)
min(rtp.homGP$RunTrue_1)
min(rtp.homGPCheap$RunTrue_1)
min(rtp.hetGPCheap$RunTrue_1)

length(which(na.omit(Run_hetML_meanvals.10mG0)<5))
length(which(na.omit(Run_homML_meanvals.10mG0)<5))
length(which(na.omit(Run_homGP_meanvals.10mG0)<5))
length(which(na.omit(Run_homGPCheap_meanvals.10mG0)<5))
length(which(na.omit(Run_hetGPCheap_meanvals.10mG0)<5))

length(which(na.omit(Run_hetML_meanvals.10mG0)<3))
length(which(na.omit(Run_homML_meanvals.10mG0)<3))
length(which(na.omit(Run_homGP_meanvals.10mG0)<3))
length(which(na.omit(Run_homGPCheap_meanvals.10mG0)<3))
length(which(na.omit(Run_hetGPCheap_meanvals.10mG0)<3))



################################################################################
################################################################################

#now do this for the mean predictions- truth

Run_homML_diff<- RunTrue.10m
values(Run_homML_diff)<- homML_mean - truevals.10m

Run_homGP_diff<- RunTrue.10m
values(Run_homGP_diff)<- homGP_mean - truevals.10m

Run_homGPCheap_diff<- RunTrue.10m
values(Run_homGPCheap_diff)<- homGPCheap_mean - truevals.10m

################################################################################
#HomML 

crs(Run_homML_diff)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homML_diff.10mG0<- Run_homML_diff
Run_homML_diffvals.10m<- raster::extract(Run_homML_diff,coords.10m)
ZeroInds<- which(Run_homML_diffvals.10m==0)
Run_homML_diffvals.10mG0<- Run_homML_diffvals.10m
Run_homML_diffvals.10mG0[ZeroInds]<- NA
values(Run_homML_diff.10mG0)<- Run_homML_diffvals.10mG0
Run_homML_diffvals.10mincm<- as.integer(round(100*Run_homML_diffvals.10mG0))
values(Run_homML_diff.10mG0)<- Run_homML_diffvals.10mincm #convert to cm since rtp can only deal with integer values

r.homMLdiff <- projectRaster(Run_homML_diff.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homMLdiff)="INT4S"

rtp.homMLdiff <- rasterToPolygons(r.homMLdiff)
rtp.homMLdiff$RunTrue_1<- rtp.homMLdiff$RunTrue_1/100


################################################################################
#HomGP50

crs(Run_homGPCheap_diff)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homGPCheap_diff.10mG0<- Run_homGPCheap_diff
Run_homGPCheap_diffvals.10m<- raster::extract(Run_homGPCheap_diff,coords.10m)
ZeroInds<- which(Run_homGPCheap_diffvals.10m==0)
Run_homGPCheap_diffvals.10mG0<- Run_homGPCheap_diffvals.10m
Run_homGPCheap_diffvals.10mG0[ZeroInds]<- NA
values(Run_homGPCheap_diff.10mG0)<- Run_homGPCheap_diffvals.10mG0
Run_homGPCheap_diffvals.10mincm<- as.integer(round(100*Run_homGPCheap_diffvals.10mG0))
values(Run_homGPCheap_diff.10mG0)<- Run_homGPCheap_diffvals.10mincm #convert to cm since rtp can only deal with integer values

r.homGPCheapdiff <- projectRaster(Run_homGPCheap_diff.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homGPCheapdiff)="INT4S"

rtp.homGPCheapdiff <- rasterToPolygons(r.homGPCheapdiff)
rtp.homGPCheapdiff$RunTrue_1<- rtp.homGPCheapdiff$RunTrue_1/100
################################################################################
#HomGP10

crs(Run_homGP_diff)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homGP_diff.10mG0<- Run_homGP_diff
Run_homGP_diffvals.10m<- raster::extract(Run_homGP_diff,coords.10m)
ZeroInds<- which(Run_homGP_diffvals.10m==0)
Run_homGP_diffvals.10mG0<- Run_homGP_diffvals.10m
Run_homGP_diffvals.10mG0[ZeroInds]<- NA
values(Run_homGP_diff.10mG0)<- Run_homGP_diffvals.10mG0
Run_homGP_diffvals.10mincm<- as.integer(round(100*Run_homGP_diffvals.10mG0))
values(Run_homGP_diff.10mG0)<- Run_homGP_diffvals.10mincm #convert to cm since rtp can only deal with integer values

r.homGPdiff <- projectRaster(Run_homGP_diff.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homGPdiff)="INT4S"

rtp.homGPdiff <- rasterToPolygons(r.homGPdiff)
rtp.homGPdiff$RunTrue_1<- rtp.homGPdiff$RunTrue_1/100
################################################################################

max(rtp.homMLdiff$RunTrue_1)
max(rtp.homGPdiff$RunTrue_1)
max(rtp.homGPCheapdiff$RunTrue_1)

min(rtp.homMLdiff$RunTrue_1)
min(rtp.homGPdiff$RunTrue_1)
min(rtp.homGPCheapdiff$RunTrue_1)

#par(mfrow=c(3,2))

maxmax<- round(max(c(max(rtp.homMLdiff$RunTrue_1),max(rtp.homGPdiff$RunTrue_1),max(rtp.homGPCheapdiff$RunTrue_1))),1)
minmin<- round(min(c(min(rtp.homMLdiff$RunTrue_1),min(rtp.homGPdiff$RunTrue_1),min(rtp.homGPCheapdiff$RunTrue_1))),1)
################################################################################
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomML_diffpred2011.jpeg",width = 800, height = 600)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homMLdiff, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homMLdiff$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Difference (m)", colors = terrain.colors(100*(2*maxmax)), limits= c(-maxmax,maxmax)) +
  labs(title= "Prediction from HomMR emulation-calibration - observation")+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))

dev.off()
beepr::beep()

################################################################################
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomGP50_diffpred2011.jpeg",width = 800, height = 600)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homGPCheapdiff, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homGPCheapdiff$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Difference (m)", colors = terrain.colors(100*(2*maxmax)), limits= c(-maxmax,maxmax)) +
  labs(title= "Prediction from HomGP50 emulation-calibration - observation") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()
beepr::beep()
################################################################################
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomGP10_diffpred2011.jpeg",width = 800, height = 600)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homGPdiff, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homGPdiff$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Difference (m)", colors = terrain.colors(100*(2*maxmax)), limits= c(-maxmax,maxmax)) +
  labs(title= "Prediction from HomGP10 emulation-calibration - observation") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()
beepr::beep()
################################################################################

#high residuals
################################################################################
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomML_diffpred2011g1.jpeg",width = 800, height = 600)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homMLdiff, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homMLdiff$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Difference (m)", colors = terrain.colors(100*(2*maxmax)), limits= c(1,maxmax)) +
  labs(title= "Prediction from HomMR emulation-calibration - Observation")
dev.off()
beepr::beep()

################################################################################
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomGP10_diffpred2011g1.jpeg",width = 800, height = 600)
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homGPdiff, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homGPdiff$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Difference (m)", colors = terrain.colors(100*(2*maxmax)), limits= c(1,maxmax)) +
  labs(title= "Prediction from HomGP10 emulation-calibration - Observation")
dev.off()
beepr::beep()

