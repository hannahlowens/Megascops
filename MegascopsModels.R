library(raster)
library(ENMeval)
library(biomod2)
library(rangeModelMetadata)

workingDirectory <- "~/Dropbox/UFBISeedGrant/Megascops/NewAnalysis/"
setwd(workingDirectory)

#Basic range model metadata information
rmmAdelpha <- rmmTemplate();
rmmAdelpha <- rmmAutofillPackageCitation(rmmAdelpha, c('raster', 'ENMeval', 
                                                       'biomod2'))
rmmAdelpha$authorship$names <- 'Owens, Hannah'
rmmAdelpha$authorship$contact <- 'hannah.owens@gmail.com'

# Load control files
regularList <- read.csv("Final_scripts/standardModels.csv", stringsAsFactors = F) 
buffer250List <- read.csv("Final_scripts/buffer250km.csv", stringsAsFactors = F)

# Get environmental data ----
setwd("ResampledChelsaAscii/")
envtList <- list.files(pattern = ".asc") #Gets a list of .asc files
envtList <- envtList[c(1,4,12,15)]
envtStack <- raster::stack(envtList, quick = F) #Reads in .asc files as a raster stack
crs(envtStack) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"; #Defines projection of layers
altitude <- crop(x = raster("../ETOPO1_Bed_g_geotiff.tif"), envtStack[[1]]);
altitude <- resample(altitude, envtStack[[1]]);
altitude <- mask(x = altitude, envtStack[[1]]);
envtStack <- stack(altitude, envtStack);
names(envtStack)[[1]] <- c("Chelsa_altitude");
plot(envtStack[[1]]); #Plots all the layers of the raster stack object
rmmAdelpha$data$environment$projection <- crs(envtStack);
rm(envtList, altitude)

# Occ files ----
setwd(paste0(workingDirectory, "/Final_Data/Input/"))
csvList <- list.files(pattern = ".csv");
csvNamesToMatch <- sub(".csv", replacement = "", x = csvList);
count <- 1;
pointList <- vector("list", length(csvList));
rmmList <- vector("list", length(csvList));
while (count <= length(pointList)){
  pointList[[count]] <- read.csv(csvList[[count]], header = T);
  
  #Remove duplicate points
  pointList[[count]] <- unique(round(cbind(pointList[[count]]$lon_dec, pointList[[count]]$lat_dec), digits = 10))
  pointList[[count]] <- pointList[[count]][complete.cases(pointList[[count]]),]
  
  #Remove points with incomplete environmental data
  pointList[[count]] <- cbind(pointList[[count]], extract(envtStack, pointList[[count]]))
  pointList[[count]] <- as.data.frame(pointList[[count]][complete.cases(pointList[[count]]),])
  pointList[[count]] <- cbind(rep(csvNamesToMatch[[count]], times = nrow(pointList[[count]])), pointList[[count]])
  colnames(pointList[[count]])[1:3] <- c("Species", "Longitude", "Latitude")
  
  #Write new files
  write.csv(pointList[[count]], paste0(workingDirectory,"/Final_Data/ProcessedLocalities/",csvList[[count]]), row.names = F)
  
  #Instantiating RMM for each species and filling in relevant data
  rmmList[[count]] <- rmmAdelpha;
  rmmList[[count]]$data$occurrence$taxon <- csvList[[count]]
  rmmList[[count]]$data$occurrence$dataType <- "Presence only"
  rmmList[[count]]$data$occurrence$sources <- "ASK KEITH"
  rmmList[[count]]$data$occurrence$presenceSampleSize <- nrow(pointList[[count]])
  
  rmmList[[count]]$dataPrep$errors$duplicateRemoval$rule <- "Coordinate";
  rmmList[[count]]$dataPrep$errors$questionablePointRemoval$notes <- "Expert-curated";
  rmmList[[count]]$dataPrep$geographic$geographicalOutlierRemoval$rule <- "Points outside training region removed";
  
  count <- count + 1
}
rm(csvList, rmmAdelpha)

# Save point rmm
saveRDS(rmmList, paste0(workingDirectory, "/Final_Data/AdelphaRMM"))

# Background files ----
rmmList <- readRDS(paste0(workingDirectory, "/Final_Data/AdelphaRMM"));

# Crop and save training regions
setwd(paste0(workingDirectory, "/Final_Data/Input/"))
shpList <- list.files(pattern = ".shp");
taxnames <- unlist(strsplit(shpList, "N.shp"));
taxnames <- unlist(strsplit(taxnames, ".shp"));
count <- 1
while (count <= length(shpList)){
  shapefile <- rgdal::readOGR(shpList[[count]]); #Reads in your shapefile
  shapeTrans <- spTransform(shapefile, "+proj=cea +lat_ts=0 +lon_0=-76.640625")
  if(taxnames[[count]] %in% buffer250List$species){
    shape2 <- buffer(x = shapeTrans, width = 250000, dissolve = T);
    rmmList[[count]]$data$environment$extentRule <- "Hull buffered by 250 km."
  }
  else{
    shape2 <- buffer(x = shapeTrans, width = 25000, dissolve = T)
    rmmList[[count]]$data$environment$extentRule <- "Hull buffered by 25 km."
  }
  shape2 <- spTransform(shape2, crs(envtStack))
  croppedStack <- crop(envtStack, shape2) #Crops raster stack to extent of training shapefile
  croppedStack <- raster::mask(croppedStack, shape2) #Masking the raster stack
  crs(croppedStack) <- crs(envtStack)
  plot(croppedStack[[1]]);
  dir.create(paste0(workingDirectory, "Final_Data/TrainingRegions/", taxnames[[count]]), 
             showWarnings = FALSE) #stops warnings if folder already exists
  setwd(paste0(workingDirectory, "Final_Data/TrainingRegions/",taxnames[[count]], "/"))
  writeRaster(croppedStack, filename = "_", format = "ascii", bylayer = T, 
              suffix=names(envtStack), NAFlag = "-9999", overwrite = T);
  
  #RMM Data
  rmmList[[count]] <- rmmAutofillEnvironment(rmmList[[count]], croppedStack, transfer = 0);
  
  #File hygene
  setwd(paste0(workingDirectory, "/Final_Data/Input/"))
  rm(croppedStack);
  
  count <- count + 1;
}

rm(shapefile, shapeTrans, shape2, count)

saveRDS(rmmList, paste0(workingDirectory, "/Final_Data/AdelphaRMM"))

# Make niche models ----
rmmList <- readRDS(paste0(workingDirectory, "/Final_Data/AdelphaRMM"));

#Get M list
setwd(paste0(workingDirectory, "Final_Data/TrainingRegions/"))
mList <- list.dirs()[c(-1)]

#Make objects for reporting
featureList <- vector(mode = "list", length = length(pointList));
regMultList <- vector(mode = "list", length = length(pointList));
deltAICList <- vector(mode = "list", length = length(pointList));
AUCdiffList <- vector(mode = "list", length = length(pointList));
thresh95List <- vector(mode = "list", length = length(pointList));
threshTModelList <- vector(mode = "list", length = length(pointList));
thresh95ModelList <- vector(mode = "list", length = length(pointList));
occPointsList <- vector(mode = "list", length = length(pointList));
trainingPointsList <- vector(mode = "list", length = length(pointList));
AUCtrain <- vector(mode = "list", length = length(pointList));
AUCtest <- vector(mode = "list", length = length(pointList));
AUCvariance <- vector(mode = "list", length = length(pointList));

count <- 1;
while (count <= length(mList)){
  print(paste("Now serving ", mList[[count]], "...", sep = ""))
  print(paste0("Indexing check: points are ", pointList[[count]][1,1]))
  setwd(mList[count])
  envtList <- list.files(pattern = ".asc")
  envtList <- envtList
  envt.st <- stack(envtList[-1])
  
  #Running the models
  if (nrow(pointList[[count]]) <= length(names(envt.st)) * 1.5){
    print("Insufficient points.")
    occPointsList[[count]] <- nrow(pointList[[count]])
    trainingPointsList[[count]] <- "None"
    featureList[[count]] <- "None"
    regMultList[[count]] <- "None"
    deltAICList[[count]] <- "None"
    AUCdiffList[[count]] <- "None"
    AUCtrain[[count]] <- "None"
    AUCtest[[count]] <- "None"
    AUCvariance[[count]] <- "None"
    thresh95List[[count]] <- "None"
    rmmList[[count]]$modelFit$algorithm <- "95% envelope"
  }#Skipping spp with a minimal number of points
  else if (nrow(pointList[[count]] <= 50)) {
    occPointsList[[count]] <- nrow(pointList[[count]])
    trainingPointsList[[count]] <- min(c((ncell(envt.st[[1]])-cellStats(envt.st[[1]], stat = "countNA"))/100, 
                                         occPointsList[[count]]*25));
    bgPoints <- randomPoints(envt.st[[1]], n = trainingPointsList[[count]], 
                             p = pointList[[count]][,2:3], excludep = T, prob = F, tryf=1.5, 
                             lonlatCorrection = T)
    trainingPointsList[[count]] <- nrow(bgPoints)
    test <- ENMevaluate(occ = pointList[[count]][,2:3], env = envt.st, bg.coords = bgPoints, 
                        RMvalues = seq(0.5,4,0.5), fc = c("L", "P", "Q", "LP", "PQ", "LQ", "LPQ"), 
                        method = "randomkfold", kfolds = 2, bin.output = T, clamp = F, parallel = T, numCores = 20, alg = "maxnet");
    minIsMax <- NULL
    numMinMax <- NULL
    check <- 1
    while(check <= length(test@predictions@layers)){
      minIsMax<- append(minIsMax, summary(test@predictions@layers[[check]])[1]==summary(test@predictions@layers[[check]])[5])
      numMinMax <- append(numMinMax, check)
      check <- check + 1
    }
    nonErrorIndex <- numMinMax[!minIsMax]
    dist <- pointDistance(cbind(test@results[nonErrorIndex,]$avg.diff.AUC, test@results[nonErrorIndex,]$AICc),c(0,0), lonlat=F)
    minDist <- min(dist)
    indexMin <- which(dist==minDist)[1]
    featureList[[count]] <- as.character(test@results[nonErrorIndex,][indexMin,]$features)
    regMultList[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$rm)
    deltAICList[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$delta.AICc)
    AUCdiffList[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$avg.diff.AUC)
    AUCtrain[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$train.AUC)
    AUCtest[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$avg.test.AUC)
    AUCvariance[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$var.test.AUC)
    predictionMap <- test@predictions[[which(names(test@predictions)==test@results[nonErrorIndex,][indexMin,1])]]
    rmmAutofillENMeval(rmmList[[count]], test, selectionCriteria = "Minimize AUCDiff and AICc scores.", 
                       optimalModelIndex = indexMin);
    
    #Thresholding
    test@occ.pts <- cbind(test@occ.pts, rep.int(1, nrow(test@occ.pts)))
    colnames(test@occ.pts) <- c("Long", "Lat", "Occurrence")
    test@bg.pts <- cbind(test@bg.pts, rep.int(0, nrow(test@bg.pts)))
    colnames(test@bg.pts) <-  c("Long", "Lat", "Occurrence")
    allPoints <- rbind(test@occ.pts, test@bg.pts)
    normPredictionMap <- (predictionMap - predictionMap@data@min)/(predictionMap@data@max-predictionMap@data@min)
    prediction <- extract(normPredictionMap, allPoints[,1:2])
    thresh95List[[count]] <- sort(extract(normPredictionMap, test@occ.pts[,1:2]), decreasing = T)[round(nrow(test@occ.pts)*.95,0)];
    thresh95ModelList[[count]] <- BinaryTransformation(normPredictionMap, thresh95List[[count]])
    writeRaster(thresh95ModelList[[count]], paste(workingDirectory, "Final_Data/RawModelOutputs/",
                                                  csvNamesToMatch[[count]],"_95ThreshModel.asc", sep = ""), 
                format = "ascii", overwrite = T);
  }
  else{
    occPointsList[[count]] <- nrow(pointList[[count]])
    trainingPointsList[[count]] <- min(c((ncell(envt.st[[1]])-cellStats(envt.st[[1]], stat = "countNA"))/100, 
                                         occPointsList[[count]]*25));
    bgPoints <- randomPoints(envt.st[[1]], n = trainingPointsList[[count]], 
                             p = pointList[[count]][,2:3], excludep = T, prob = F, tryf=1.5, 
                             lonlatCorrection = T)
    test <- ENMevaluate(occ = pointList[[count]][,2:3], env = envt.st, bg.coords = bgPoints, 
                        RMvalues = seq(0.5,4,0.5), fc = c("L", "P", "Q", "LP", "PQ", "LQ", "LPQ"), 
                        method = "block", bin.output = T, clamp = F, parallel = T, numCores = 20, alg = "maxnet");
    minIsMax <- NULL
    numMinMax <- NULL
    check <- 1
    while(check <= length(test@predictions@layers)){
      minIsMax<- append(minIsMax, summary(test@predictions@layers[[check]])[1]==summary(test@predictions@layers[[check]])[5])
      numMinMax <- append(numMinMax, check)
      check <- check + 1
    }
    nonErrorIndex <- numMinMax[!minIsMax]
    dist <- pointDistance(cbind(test@results[nonErrorIndex,]$avg.diff.AUC, test@results[nonErrorIndex,]$AICc),c(0,0), lonlat=F)
    minDist <- min(dist)
    indexMin <- which(dist==minDist)[1]
    featureList[[count]] <- as.character(test@results[nonErrorIndex,][indexMin,]$features)
    regMultList[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$rm)
    deltAICList[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$delta.AICc)
    AUCdiffList[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$avg.diff.AUC)
    AUCtrain[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$train.AUC)
    AUCtest[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$avg.test.AUC)
    AUCvariance[[count]] <- as.numeric(test@results[nonErrorIndex,][indexMin,]$var.test.AUC)
    predictionMap <- test@predictions[[which(names(test@predictions)==test@results[nonErrorIndex,][indexMin,1])]]
    rmmAutofillENMeval(rmmList[[count]], test, selectionCriteria = "Minimize AUCDiff and AICc scores.", 
                       optimalModelIndex = indexMin);
    
    #Thresholding
    test@occ.pts <- cbind(test@occ.pts, rep.int(1, nrow(test@occ.pts)))
    colnames(test@occ.pts) <- c("Long", "Lat", "Occurrence")
    test@bg.pts <- cbind(test@bg.pts, rep.int(0, nrow(test@bg.pts)))
    colnames(test@bg.pts) <-  c("Long", "Lat", "Occurrence")
    allPoints <- rbind(test@occ.pts, test@bg.pts)
    normPredictionMap <- (predictionMap - predictionMap@data@min)/(predictionMap@data@max-predictionMap@data@min)
    prediction <- extract(normPredictionMap, allPoints[,1:2])
    thresh95List[[count]] <- sort(extract(normPredictionMap, test@occ.pts[,1:2]), decreasing = T)[round(nrow(test@occ.pts)*.95,0)];
    thresh95ModelList[[count]] <- BinaryTransformation(normPredictionMap, thresh95List[[count]])
    writeRaster(thresh95ModelList[[count]], paste(workingDirectory, "Final_Data/RawModelOutputs/", 
                                                  csvNamesToMatch[[count]],"_95ThreshModel.asc", sep = ""), format = "ascii", overwrite = T);
  }
  
  setwd(paste0(workingDirectory, "Final_Data/TrainingRegions/"))
  count <- count + 1;
}

saveRDS(rmmList, paste0(workingDirectory, "/Final_Data/AdelphaRMM"))


reportTable <- as.data.frame(cbind(unlist(occPointsList), unlist(trainingPointsList), unlist(featureList),
                                   unlist(regMultList),unlist(deltAICList), unlist(AUCdiffList), unlist(AUCtrain), 
                                   unlist(AUCtest), unlist(AUCvariance), unlist(thresh95List)));
colnames(reportTable) <- c("Number of Occurrences","Number of Pseudoabsences","Features", "Regularization Multiplier", 
                           "Delta AICc", "Mean AUC Difference", "AUC Train", "Mean AUC Test", "AUC Variance", "Threshold95")
rownames(reportTable) <- taxnames
write.csv(as.matrix(reportTable), paste(workingDirectory, "Final_Data/RawModelOutputs/ResultsTable.csv", sep = ""));
