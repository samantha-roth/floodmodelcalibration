flood_extent_run<-function(i){
  library(raster)
  outputs <- data.frame()
  river <- sample_river
  par <- sample_par
  
  #river[3, 7] <- 574600*0.3048^3 # FEMA 100-yr upstream discharge in cms
  river[3, 7] <- 13110.70 # 2011 flood upstream discharge in cms
  
  river[3, 4] <- mat[i,"n_ch"] # upstream channel roughness
  river[7, 4] <- mat[i,"n_ch"] # downstream channel roughness
  ##STUFF I AM ADDING##
  river[3:7, 3] <- as.character(as.numeric(river[3:7, 3])*mat[i,"rwe"]) # river width * error
  #river[c(3,7), 5] <- river[c(3,7), 5]+ mat[i,"ree"] # river bed elevation + error 
  river[3, 5] <- as.character(as.numeric(river[3, 5])+ mat[i,"ree"]) # river bed elevation + error 
  river[7, 5] <- as.character(as.numeric(river[7, 5])+ mat[i,"ree"]) # river bed elevation + error 
  #####################
  par[6,2] <- mat[i,"n_fp"] # floodplain roughness
  par[5,2] <- DEM_file
  
  par[1, 2] <- paste0("run",i)
  
  setwd(paste0(wd,"/LISFLOOD"))
  write.table(
    river,
    "./Selinsgrove.river",
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t",
    na = ""
  )
  write.table(
    par,
    "./Selinsgrove.par",
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t",
    na = ""
  )
  
  print(paste("Run", i))
  # Run LISFLOOD-FP on Desktop
  #system("./lisflood -v Selinsgrove.par")
  # Run LISFLOOD-FP on cluster
  system("./lisflood.exe Selinsgrove.par")
  file.rename(paste0("./results/run", i, ".max"), "./results/max.asc")
  
  #flood extent
  model_hazard <- raster("./results/max.asc", format = "ascii")
  
  # saving results: model flood extend, hazard and risk
  setwd(wd)
  
  outputs[1, 1] <- i
  outputs[1, 2] <- mat[i,"n_ch"]
  outputs[1, 3] <- mat[i,"n_fp"]
  outputs[1, 4] <- mat[i,"rwe"]
  outputs[1, 5] <- mat[i,"ree"]
  
  
  # Save results

if(prior==1){

if(set=="first"){

  	writeRaster(model_hazard,
    	paste0("./Outputs10m/prior1/Extent/Run_", i, ".asc"),
   	format = "ascii",
    	overwrite = T)
  
  	colnames(outputs) <-c("run no.","n_ch","n_fp","rwe","ree")
  	write.csv(outputs, paste0("./Outputs10m/prior1/Run_", i, ".csv"))

   }
   

   if(set=="second"){
  	writeRaster(model_hazard,
    	paste0("./Outputs10m/prior1/second/Extent/Run_", i, ".asc"),
    	format = "ascii",
    	overwrite = T)
  
  	colnames(outputs) <-c("run no.","n_ch","n_fp","rwe","ree")
  	write.csv(outputs, paste0("./Outputs10m/prior1/second/Run_", i, ".csv"))

   }

   if(set=="third"){
  	writeRaster(model_hazard,
    	paste0("./Outputs10m/prior1/third/Extent/Run_", i, ".asc"),
    	format = "ascii",
    	overwrite = T)
  
  	colnames(outputs) <-c("run no.","n_ch","n_fp","rwe","ree")
  	write.csv(outputs, paste0("./Outputs10m/prior1/third/Run_", i, ".csv"))

   }

   if(set=="fourth"){
  	writeRaster(model_hazard,
    	paste0("./Outputs10m/prior1/fourth/Extent/Run_", i, ".asc"),
    	format = "ascii",
    	overwrite = T)
  
  	colnames(outputs) <-c("run no.","n_ch","n_fp","rwe","ree")
  	write.csv(outputs, paste0("./Outputs10m/prior1/fourth/Run_", i, ".csv"))

   }
  
  if(set=="fifth"){
    writeRaster(model_hazard,
                paste0("./Outputs10m/prior1/fifth/Extent/Run_", i, ".asc"),
                format = "ascii",
                overwrite = T)
    
    colnames(outputs) <-c("run no.","n_ch","n_fp","rwe","ree")
    write.csv(outputs, paste0("./Outputs10m/prior1/fifth/Run_", i, ".csv"))
    
  }
  
  if(set=="sixth"){
    writeRaster(model_hazard,
                paste0("./Outputs10m/prior1/sixth/Extent/Run_", i, ".asc"),
                format = "ascii",
                overwrite = T)
    
    colnames(outputs) <-c("run no.","n_ch","n_fp","rwe","ree")
    write.csv(outputs, paste0("./Outputs10m/prior1/sixth/Run_", i, ".csv"))
    
  }
  
  if(set=="seventh"){
    writeRaster(model_hazard,
                paste0("./Outputs10m/prior1/seventh/Extent/Run_", i, ".asc"),
                format = "ascii",
                overwrite = T)
    
    colnames(outputs) <-c("run no.","n_ch","n_fp","rwe","ree")
    write.csv(outputs, paste0("./Outputs10m/prior1/seventh/Run_", i, ".csv"))
    
  }
 

}

if(prior==2){

  writeRaster(
    model_hazard,
    paste0("./Outputs10m/prior2/Extent/Run_", i, ".asc"),
    format = "ascii",
    overwrite = T
  )
  
  
  colnames(outputs) <-
    c("run no.","n_ch","n_fp","rwe","ree")
  write.csv(outputs, paste0("./Outputs10m/prior2/Run_", i, ".csv"))
}




}
