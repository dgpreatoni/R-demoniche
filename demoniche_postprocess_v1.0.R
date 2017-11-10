################################################################################
# demoniche output postprocessor
# 
# version 1.0
# created: prea 20170504
# updated: elia 20170505 - created exerciser code
#          prea 20170522 - some cleanup and reordering
#          elia 20170530 - add values in "demoniche_setup"
#          elia 20170531 - add values in "demoniche_setup"
#          elia 20170601 - create matrix MastYear and PoorYear
#          prea 20170607 - modified to run parallel on eira
#          prea 20170710 - reworked as a stand-alone script
################################################################################

# clean up R environment
rm(list=ls())

# invoke the packages we need
library(readxl) # read Excel files
library(rgdal) # read/write shapefiles
library(raster) # read/write/manage rasters
#library(demoniche) # also loads sp
library(car) # recode
library(ggplot2)
library(ggmap)
require(gtools) #quantcut
library(RColorBrewer)

################################################################################
# user-defined variables here!

# define some variables to store directory names
rootDir <- '/data/RData/demoniche_test'

# .rda file created by demoniche_setup (assumed in <rootDir>/output)
modelFile <- 'test_1500_v2.32norm.rda'

# the data directory defined by 'foldername' when demoniche_model() has been run (assumed in <rootDir>/outout)
scenarioDir <- 'test_1500_v2.32norm'

# end of user-defined variables block
################################################################################


# build subdirectories names from rootDir
baseDir <- paste(rootDir, 'BASE', sep='/')
dataDir <- paste(rootDir, 'data', sep='/')
outputDir <- paste(rootDir, 'output', sep='/')
scenarioDir <- paste(outputDir, scenarioDir, sep='/')

# load utilities for meters<->degrees conversions
source(paste0(rootDir,'/scripts/meters-to-degrees.R'))

# load model data: be aware that elements named "dispersal_probabilities" and "dist_latlong" can be _HUGE_ (and we don't need them here)
## function to load R datafile and return read objects as to place them into a declared variable
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
model <- loadRData(paste(outputDir, modelFile, sep='/'))
# get rid of 'heavy' elements
model$dispersal_probabilities <- NULL
model$dist_latlong <- NULL
gc() # call explicitly the garbage collector to free up RAM!
# all the parametrizations used in modeling are now available in 'model' object, except dispersal probabilities and site distance matrix

# read in the "reference raster"
refRast.ids <- raster(paste0(scenarioDir,'_reference_raster'))
refRast <- refRast.ids
values(refRast) <- 0

#### reload output data, try making out maps
setwd(scenarioDir)

#### SANITY CHECKS HERE

# write out reference raster
refRast.ids.sp <- rasterToPolygons(refRast.ids)
names(refRast.ids.sp) <- 'cellID'
writeOGR(refRast.ids.sp, dsn=getwd(), layer="Reference_raster", driver="ESRI Shapefile", overwrite=TRUE)

# write out initial populations

Orig_Populations.sp <- SpatialPointsDataFrame(coords=model$Orig_Populations[,c('XCOORD', 'YCOORD')], data=model$Orig_Populations[,c('PatchID', 'area_population')], proj4string=CRS(projection(refRast)))
names(Orig_Populations.sp) <- c("PatchID", "area_pop")
Orig_Populations.sp$dens <- model$density_individuals
Orig_Populations.sp$Ncalc <- Orig_Populations.sp$dens * Orig_Populations.sp$area_pop
writeOGR(Orig_Populations.sp, dsn=getwd(), layer="Orig_Populations", driver="ESRI Shapefile", overwrite=TRUE)

#### END OF SANITY CHECKS

# grab all Projection repetition files, make out a list where each element is a replicate
projFiles <- list.files(pattern='Projection_rep')      
Projections <- list()
for(f in projFiles) {
  load(f)
  Projections[[substr(f,1,nchar(f)-4)]] <- Projection
}
rm(f)
rm(Projection)

#### Notes:
# An element of 'Projections' is a plain 4-dimensional R named array. Array 'payload' is number of individuals (not sure...)
# - Dimension 1 is _time_, i.e. the number of year steps defined as an integer in no_yrs parameter in demoniche_setup(). Dimension names are 'timesliceyear_<n>'.
# - Dimension 2 is _stage_, i.e. the number of stages defined with the 'stages' parameter in demoniche_setup(). Dimension names are identical to 'stages'.
# - Dimension 3 is _patch ID_, identical to those in the gridID column of the Nichemap dataframe supplied in demoniche_setup().
# - Dimension 4 is _scenario matrix_, i.e one value for each scenario matrix defined in the matrices parameter in demoniche_setup().
# - Summing up: aPrj[<timeslice index>, <stage index>, <patch ID>, <scenario>]. Index 3 (patchID) should be always left blank to extract all patches.
# - To have a 'flat' matrix with number of individuals, a 'Projection' must be sliced and multiplied by the 'sumWeight' parameter, and columns must be summed up, like this: colSums(aPrj[<timeSlice>,,,<scenario>] * sumWeight)

#### define, if need be, accessor functions here
# returns how many years a projection spans
lengthYears <- function(aPrj) {
  # also works for lists of replicas...
  if(inherits(aPrj, 'list')) { # then is a list, get data from first element
    return(as.numeric(lapply(aPrj, function(x) lengthYears(x)))) # iteration is human, recursion is divine...
  } else { # not a list
    return(dim(aPrj)[1])
  }
}
# returns valid ordered names for dimension 4 (scenarios)
scenarioNames <- function(aPrj) {
  return(dimnames(aPrj)[[4]])
}
# returns a number of individuals dataframe at a given time slice
getIndividualsAtYear <- function(aPrj, aYear, sumWeight, aScenario=1) {
  N <- colSums(aPrj[aYear,,,aScenario] * sumWeight)
  patchID <- as.integer(dimnames(aPrj)[[3]])
  return(data.frame(patchID=patchID,N=N))
}

# read in and reshape data
projectionsMaps <- list()
projectionsPopulations <- data.frame()
# cycle through replicates
for(p in names(Projections)) {
  numYears <- lengthYears(Projections[[p]])
  Scenarios <- scenarioNames(Projections[[p]])
  scenariosMaps <- list()
  # create a stack for each scenario
  for(s in Scenarios) {
    scenarioData <- data.frame()
    mapStack <- stack()
    for(y in 1:numYears) {
      df <- getIndividualsAtYear(Projections[[p]], y, model$sumweight, s)
      #write.csv(df, file=paste0('individuals-scenario_', s, '_rep', p, '_year', sprintf("%04d", y), '.csv'))
      mapY <- refRast
      mapY[df$patchID] <- df$N
      names(mapY) <- paste('Year', sprintf("%04d", y), sep='_')
      mapStack <- addLayer(mapStack, mapY)
      names(df) <- c('patchID', paste0('Y', sprintf("%04d", y)))
      if(y==1){
        scenarioData <- df
      } else {
        scenarioData <- merge(scenarioData, df, by='patchID')
      }
      rm(df, mapY)
    }
    scenariosMaps[[s]] <- mapStack
    scenarioData$Scenario <- s 
    scenarioData$Replica <- p
    projectionsPopulations <- rbind(projectionsPopulations, scenarioData)
    # this dumps a 'raw' stack for the current scenario
    #writeRaster(mapStack, filename=paste0(substr(modelFile, 1, nchar(modelFile)-4), '-Scenario_', s, '.img'), format="HFA", overwrite=TRUE)
    #writeRaster(mapStack, filename=paste0(p, '-Scenario_', s, '.img'), format="HFA", overwrite=TRUE)
    rm(mapStack, scenarioData)
  }
  projectionsMaps[[p]] <- scenariosMaps
  rm(scenariosMaps)
}
# projectionsMaps now contains all the rasters.
# projectionsMaps has the following structure:
# projectionMaps[[<replica name>]][[<scenario name]] is a raster::stack, one layer for each year named "Year_<nnnn>"
# projectionsPopulations now contains all the numerical estimates per year and cell, tagged by scenario name and replica name

# now, aggregate replicates, keeping scenarios separated
averagedMaps <- list()
densityMaps <- list()
numYears <- lengthYears(Projections[[1]])
for(s in Scenarios) {
  cat("Scenario:", s, "\n")
  scenarioStack <- stack()
  for(y in 1:numYears) {
    cat("\tyear:", y, '\n')
    replicasStack <- stack()
    for(p in names(Projections)[1:3]) {
      cat("\t\treplica:", p, '\n')
      replicasStack <- addLayer(replicasStack, projectionsMaps[[p]][[s]][[paste('Year', sprintf("%04d", y), sep='_')]])
    }
    yearRaster <- mean(replicasStack)
    names(yearRaster) <- paste('Year', sprintf("%04d", y), sep='_')
    rm(replicasStack)
    scenarioStack <- addLayer(scenarioStack, yearRaster)
  }
  # clean up map
  averagedMaps[[s]] <- clamp(scenarioStack, lower=0)
  # calculate densities
  densityMaps[[s]] <- averagedMaps[[s]] / (getAreaMeters(averagedMaps[[s]]) / 10000) # individuals per hectare 
  # write rasters for the current scenario
  writeRaster(averagedMaps[[s]], filename=paste0('Distribution-individuals-Scenario_', s, '.img'), format="HFA", overwrite=TRUE)
  writeRaster(densityMaps[[s]], filename=paste0('Distribution-density-Scenario_', s, '.img'), format="HFA", overwrite=TRUE)
  # write separated rasters, one per year
  for(l in 1:nlayers(densityMaps[[s]])) {
    writeRaster(densityMaps[[s]][[l]], filename=paste0('Distribution-density-Scenario_', s, '-Year_', sprintf("%04d", l), '.tif'), format="GTiff", overwrite=TRUE)
    writeRaster(densityMaps[[s]][[l]], filename=paste0('Distribution-individuals-Scenario_', s, '-Year_', sprintf("%04d", l), '.tif'), format="GTiff", overwrite=TRUE)
  }
}
# now averagedMaps contains one stack per scenario with raw number of animals
# and densityMaps contains one stack per scenario with population density in individuals per hectare

# try preparing stand-alone images
densityBreaks <- c( 0.01,   0.1,   0.5,   1,   2,   3,   5,    20, 100) 
densityLabels <- c('<0.01', '0.1', '0.5', '1', '2', '3', '5', '>20')
densityAlpha <-  c( 0.1,  0.2,   0.5, 0.8, 0.8, 0.8, 0.8,   0.8)
densityPalette <- c('#FBFBFB', brewer.pal(9, 'RdPu'))
#plot(densityMaps[[1]][[1]], breaks=densityBreaks, col=densityPalette)

# use ggmap, output separate PNGs per year
le <- extent(refRast)
locationExtent <- c(le[1], le[3], le[2], le[4])
background <- get_map(location=locationExtent)

#### superimpose a raster with geom_raster
makeYearMap <- function(filename, aBackground, rast, breaks, labels, title) {
  # raster must be turned into a dataframe with x, y, z.
  rast.df <- as.data.frame(rast, xy=TRUE)
  names(rast.df)[3] <- 'z'
  # pre-cut data on fixed breaks
  rast.df.cut <- data.frame(x=rast.df$x, y=rast.df$y, z=cut(rast.df$z, breaks=breaks, labels=labels, include.lowest=FALSE))
  #rast.df.cut$alpha <- as.numeric(as.character(recode(rast.df.cut$z, "'0'=0; '<0.01'=0.0; '0.1'=0.0; '0.5'=0.0; '1'=0.1; '2'=0.1; '3'=0.1; '5'=0.1; '>20'=0.1")))
  rast.df.cut <- rast.df.cut[complete.cases(rast.df.cut),]
  png(filename=filename, width=200, height=190, units='mm', res=300)
  print(ggmap(aBackground, extent='device') + geom_raster(data=rast.df.cut, aes(x=x, y=y, fill=z), alpha=0.7, interpolate=TRUE) + scale_fill_brewer(palette='RdPu', direction=1, name=title) + coord_cartesian())
  dev.off()
}
# plot (only scenario 1)
for(l in 1:nlayers(densityMaps[[1]])) {
  makeYearMap(filename=paste0('density_year-', sprintf("%04d", l), '.png'), aBackground=background, rast=densityMaps[[1]][[l]], breaks=densityBreaks, labels=densityLabels, 'ind/ha')}

## make a video
system('ffmpeg -y -r 2 -i density_year-%04d.png -b 20M output.mp4')

# plot, no control on color
#ggmap(background) + inset_raster(as.raster(rast), xmin=extent(rast)[1], xmax=extent(rast)[2], ymin=extent(rast)[3], ymax=extent(rast)[4]) + coord_map()
# plot, full control on color and whatever else
#ggmap(background) + geom_raster(data=rast.df, aes(x=x, y=y, fill=z), interpolate=TRUE) + scale_fill_distiller(palette='RdPu', direction=1) + coord_cartesian()

# now, sum up number of animals per scenario and average replicates per year
for(s in Scenarios) {
  popData <- split(projectionsPopulations[projectionsPopulations$Scenario==s,], projectionsPopulations$Replica)
  popData <- lapply(popData, function(x) colSums(x[,2:(ncol(x)-2)]))
  popData <- data.frame(t(do.call('rbind', popData)))
  popDataRaw <- data.frame()
  for(r in 1:ncol(popData)) {
    popDataRaw <- rbind(popDataRaw, data.frame(Year=row.names(popData), N=popData[,r], replica=r))
  }
  popDataRaw$Year <- as.numeric(substr(popDataRaw$Year,2,5))
  png(filename=paste0('population_trend_Scenario-', s, '.png'), width=200, height=190, units='mm', res=300)
  print(ggplot(popDataRaw[popDataRaw$replica %in% c(1,2,3),], aes(x=Year, y=N)) + geom_smooth())
  dev.off()
}
