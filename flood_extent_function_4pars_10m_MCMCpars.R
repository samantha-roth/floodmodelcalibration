flood_extent_run<-function(i){
  library(raster)
  outputs <- data.frame()
  river <- sample_river
  par <- sample_par
  
  #river[3, 7] <- 574600*0.3048^3 # FEMA 100-yr upstream discharge in cms
  #river[3, 7] <- 13110.70 # 2011 flood upstream discharge in cms
  river[3, 7] <- 12091.29 # 2004 flood upstream discharge in cms
  
  river[3, 4] <- mat[i,1] # upstream channel roughness
  river[7, 4] <- mat[i,1] # downstream channel roughness
  ##STUFF I AM ADDING##
  river[3:7, 3] <- as.character(as.numeric(river[3:7, 3])*mat[i,3]) # river width * error
  #river[c(3,7), 5] <- river[c(3,7), 5]+ mat[i,"ree"] # river bed elevation + error 
  river[3, 5] <- as.character(as.numeric(river[3, 5])+ mat[i,4]) # river bed elevation + error 
  river[7, 5] <- as.character(as.numeric(river[7, 5])+ mat[i,4]) # river bed elevation + error 
  #####################
  par[6,2] <- mat[i,2] # floodplain roughness
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
  outputs[1, 2] <- mat[i,1]
  outputs[1, 3] <- mat[i,2]
  outputs[1, 4] <- mat[i,3]
  outputs[1, 5] <- mat[i,4]
  
  
  # Save results

if(emulator=="hetML"){

  writeRaster(
    model_hazard,
    paste0("./Outputs10m/hetML/Extent/Run_", i, ".asc"),
    format = "ascii",
    overwrite = T
  )
  
  
  colnames(outputs) <-
    c("run no.","n_ch","n_fp","rwe","ree")
  write.csv(outputs, paste0("./Outputs10m/hetML/Run_", i, ".csv"))
}

if(emulator=="homML"){

  writeRaster(
    model_hazard,
    paste0("./Outputs10m/homML/Extent/Run_", i, ".asc"),
    format = "ascii",
    overwrite = T
  )
  
  
  colnames(outputs) <-
    c("run no.","n_ch","n_fp","rwe","ree")
  write.csv(outputs, paste0("./Outputs10m/homML/Run_", i, ".csv"))
}

if(emulator=="hetGP"){

  writeRaster(
    model_hazard,
    paste0("./Outputs10m/hetGP/Extent/Run_", i, ".asc"),
    format = "ascii",
    overwrite = T
  )
  
  
  colnames(outputs) <-
    c("run no.","n_ch","n_fp","rwe","ree")
  write.csv(outputs, paste0("./Outputs10m/hetGP/Run_", i, ".csv"))
}

if(emulator=="homGP"){

  writeRaster(
    model_hazard,
    paste0("./Outputs10m/homGP/Extent/Run_", i, ".asc"),
    format = "ascii",
    overwrite = T
  )
  
  
  colnames(outputs) <-
    c("run no.","n_ch","n_fp","rwe","ree")
  write.csv(outputs, paste0("./Outputs10m/homGP/Run_", i, ".csv"))
}


if(emulator=="hetGPCheap"){

  writeRaster(
    model_hazard,
    paste0("./Outputs10m/hetGPCheap/Extent/Run_", i, ".asc"),
    format = "ascii",
    overwrite = T
  )
  
  
  colnames(outputs) <-
    c("run no.","n_ch","n_fp","rwe","ree")
  write.csv(outputs, paste0("./Outputs10m/hetGPCheap/Run_", i, ".csv"))
}

if(emulator=="homGPCheap"){

  writeRaster(
    model_hazard,
    paste0("./Outputs10m/homGPCheap/Extent/Run_", i, ".asc"),
    format = "ascii",
    overwrite = T
  )
  
  
  colnames(outputs) <-
    c("run no.","n_ch","n_fp","rwe","ree")
  write.csv(outputs, paste0("./Outputs10m/homGPCheap/Run_", i, ".csv"))
}




}
