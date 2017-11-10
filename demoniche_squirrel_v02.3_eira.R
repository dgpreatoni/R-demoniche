################################################################################
# demoniche test
# 
# version 2.3
# created: prea 20170504
# updated: elia 20170505 - created exerciser code
#          prea 20170522 - some cleanup and reordering
#          elia 20170530 - add values in "demoniche_setup"
#          elia 20170531 - add values in "demoniche_setup"
#          elia 20170601 - create matrix MastYear and PoorYear
#          prea 20170607 - modified to run parallel on eira
#          elia 20170619 - create Weighted Mean matrix
#          elia 20170626 - reviewing values
#          prea 20170629 - modified to run on eira
#          elia 20170704 - now 'resanmples' the landcover raster summing up carrying capacities
#          prea 20170706 - modified to use resample functions from the 'velox' package
#                        - checked resampling technique
################################################################################

# clean up R environment
rm(list=ls())

# invoke the packages we need
library(readxl) # read Excel files
library(rgdal) # read/write shapefiles
library(raster) # read/write/manage rasters
library(velox) # speeds up raster stuff
library(demoniche) # also loads sp


# define some variables to store directory names
rootDir <- '/srv/data/demoniche_test'

# the spatial resolution at which we want to work
resolution <- 2000 # meters

# raster data with land cover map (must be in BASE directory)
landCoverFile <- 'CLC2012_4326.tif'
# excel spreadsheet with lookup table to reclassify land cover map (must be in data directory)
landCoverReclassFile <- 'clc_legend_reclass.xls'
# shapefile containing species point locations (must be in BASE directory)
sightingsFile <- 'gray_squirrel_sightings_4326.shp'

# build subdirectories names from rootDir
baseDir <- paste(rootDir, 'BASE', sep='/')
dataDir <- paste(rootDir, 'data', sep='/')
outputDir <- paste(rootDir, 'output', sep='/')
# folder name to be used to store model outputs (will be in outoputDir)
scenarioDir <- "test_2000"

# change the working directory to outputdir
setwd(outputDir)

# calculate cell size in hectares
cellSizeHa <- resolution^2/10000

# start a multicore cluster, use on eira only
beginCluster(7)

#### process land cover raster to obtain habitat suitability data frame

## read in raw landcover data and reclassify to have individuals
# read in raster
landCover <- raster(paste(baseDir, landCoverFile, sep='/'))
# read in reclass table
landCoverReclass <- read_excel(paste(dataDir, landCoverReclassFile, sep='/'))
# reformat reclass table as needed by raster::reclassify()
rclTable <- as.matrix(landCoverReclass[, c('GRID_CODE', 'GREY_K')])
# reclassify
landCover.rcl <- reclassify(landCover, rclTable)
# we have densities here, multiply them by cell area so we have individuals
# calculate cell area (of top-left) cell
cellArea <- pointDistance(coordinates(landCover.rcl)[1,], coordinates(landCover.rcl)[2,], lonlat=TRUE)^2 / 10000 # hectares!
# calculate number of individuals
landCover.rcl <- landCover.rcl * cellArea
# reproject if need be
if(!compareCRS(landCover.rcl, CRS("+init=epsg:4326"))) {
  landCover.rcl <- projectRaster(landCover.rcl, CRS("+init=epsg:4326"))
}

## resample siutability data at a coarser resolution
# the resample command needs a "reference" raster, approximate resolution in degrees 
refRast <- raster(crs=CRS("+init=epsg:4326"), ext=extent(landCover.rcl), resolution=resolution*0.00011111/12.35, vals=0)
# create a copy of the reference raster refRast having cell IDs as values
refRast.ids <- refRast
values(refRast.ids) <- 1:ncell(refRast.ids)
# convert it to polygons, since raster::extent needs polygons
refPoly <- rasterToPolygons(refRast.ids)

#'resample' carrying capacities raster summing up carrying capacities
# resampling is carried out using the velox package: we first have to convert the raster to be extracted into a velox raster
landCover.rcl.vx <- velox(landCover.rcl)
kappa <- landCover.rcl.vx$extract(sp=refPoly, fun=sum) # FASTER!!
kappa[is.na(kappa)] <- 0 # pad unsuitable areas having NA with 0

# rescale as a 0-1 "suitability" index as needed by demoniche
suitability <- kappa/max(kappa)

# prepare the niche dataframe for demoniche
nicheCoordinates <- coordinates(refRast)
nicheDataFrame <- data.frame(gridID=values(refRast.ids), X=nicheCoordinates[,'x'], Y=nicheCoordinates[,'y'], Years00=suitability, Years10=suitability, Years20=suitability, Years30=suitability)
rm(suitability, nicheCoordinates)

# clean up
rm(landCoverReclass, rclTable, landCover, landCoverFile, landCoverReclassFile)


#### process distribution data to obtain distribution data frame

# read in point locations
sightings <- readOGR(paste(baseDir, sightingsFile, sep = '/'))
if(!compareCRS(sightings, CRS("+init=epsg:4326"))) {
  sightings <- spTransform(sightings, CRS("+init=epsg:4326"))
}
# convert point locations to raster data using the reference raster created before (for this reason we calculate niche first and distribution second, since calculating niche implies producing the reference raster)
sightingsDataFrame <- extract(refRast, sightings, cellnumbers=TRUE, df=TRUE)
# aggregate counts by grid cell
sightingsDataFrame <- as.data.frame(xtabs(~cells, sightingsDataFrame))
# get cells coordinates from reference raster
sightingsCoordinates <- as.data.frame(coordinates(refRast))
sightingsCoordinates$cells <- row.names(sightingsCoordinates)
# merge coordinates with aggregated signtihgs
sightingsDataFrame <- merge(sightingsDataFrame, sightingsCoordinates, by='cells', all.x=TRUE)
# reorder the resulting dataframe as demoniche needs
names(sightingsDataFrame) <- c('patchID', 'area', 'X', 'Y')
sightingsDataFrame$patchID <- as.integer(as.character(sightingsDataFrame$patchID))
sightingsDataFrame <- sightingsDataFrame[c(1, 3, 4, 2)]
# use cells ID to read kappa values
maxIndividuals <- kappa[sightingsDataFrame$patchID]
# clean up
rm(sightingsFile, sightingsCoordinates)

#### create all the non-spatial data

# transition matrix, see for example /demoniche_test/data/Model parameters v5/Dynamics -> tab "parameters" row 3, 5, 6 col 3
matrixTrans1 <- matrix(c(0.3, 0.285, 1.302, 0.6), nrow=2, ncol=2, byrow=FALSE, dimnames = list(c("J", "A"), c("J", "A")))
matrixTrans2 <- matrix(c(0.25, 0.285, 0.5, 0.55), nrow=2, ncol=2, byrow=FALSE, dimnames=list(c("J", "A"), c("J", "A")))
matrixTrans3 <- matrix(c(0.4, 0.285, 2.45, 0.65))
matrixTrans4 <- matrix(c(0.3125, 0.35, 1.3885, 0.6))
# the matrix need to be 'unfolded' into a single column
matrixTransNorm <- as.vector(matrixTrans1)
matrixTransPoor <- as.vector(matrixTrans2)
matrixTransMast <- as.vector(matrixTrans3)
matrixTransWM <- as.vector(matrixTrans4)   # Weighted Mean matrix

# assemble transition matrixes
#matrixTrans <- matrix(c(matrixTransNorm, matrixTransPoor, matrixTransMast), nrow=4, ncol=3, dimnames=list(c("JJ", "JA", "AJ", "AA"), c("Norm", "Poor", "Mast")))
matrixTrans <- matrix((matrixTransWM), nrow=4, ncol=1, dimnames=list(c("JJ", "JA", "AJ",       "AA"), ("GreySquirrel")))

# define labels for stages
growthStages <-  c("J", "A")

# initial proportion of stages, assuming 1:1 sex ration and an average number of 5  newborn per litter
propInitial <-  c(0.39, 0.61)

# count all life stages when sizing up populations
sumWeight <- c(1,1)

# proportion of effect of carrying capacities on life stages
kWeight <- c(0,0.95)

# number of years in a time frame
noYrs <- 10

# gaussian noise 'variation matrix' (matrix_var)
matrixVar <- matrix(c(0.01, 0.05),
                    nrow = nrow(matrixTrans1),
                    ncol = 1,
                    dimnames = list(NULL, "sd"))

# short distance dispersal kernel
sDistanceDisp <- 0.05

# long ditance dispersal parameters
lDistanceDisp <- 5 # kilometers

# distance dispersion max
dispConst <- c(0.7, 0.7, 0.1, 0.045) #5000m = 0.045°

# how to hamdle scenario matrixes (probability to draw from)
probScenario <- c(0.5, 0.25, 0.25) 

# effects on demography
transAffDemogr <- "all"

# calculate initial population densities
# original densities are in sightingsDataframe$area
densityIndividuals <- sightingsDataFrame$area # / cellSizeHa
# since this is a raster-based model and not a patch-based one, replace SightingsDataframe$area with grid cell area
sightingsDataFrame$area <- cellSizeHa


#### prepare data with the setup function
demoniche_setup(modelname="GreySquirrel", 
                Populations=sightingsDataFrame,
                stages=growthStages,
                Nichemap=nicheDataFrame,
                matrices=matrixTrans,
                matrices_var=matrixVar,
                prob_scenario=probScenario,
                proportion_initial=propInitial,
                density_individuals=densityIndividuals,
                transition_affected_niche= "all",
                transition_affected_env= "all",
                transition_affected_demogr= transAffDemogr,
                env_stochas_type= "normal",
                noise=1, # not supported
                fraction_SDD=sDistanceDisp, 
                fraction_LDD=lDistanceDisp, 
                dispersal_constants=dispConst,
                no_yrs=noYrs,
                Ktype="ceiling",
                K=maxIndividuals,
                Kweight= kWeight,
                
                sumweight = sumWeight)
# run a basic simulation
test0 <- demoniche_model(modelname = "GreySquirrel",
                         Niche = TRUE,
                         Dispersal = TRUE,
                         repetitions = 5,
                         foldername = scenarioDir)


# use on eira only, close up cluster
endCluster()


#### stop here, what follows will go in a separate  post-processing script


#### reload outout data, try making out maps
setwd(paste(outputDir, scenarioDir, sep='/'))

projFiles <- list.files(pattern='Projection_rep')      
Projections <- list()
for(f in projFiles) {
  load(f)
  Projections[[f]] <- Projection
}
rm(Projection)

mapList <- list()
for(f in projFiles) {
  numYears <- dim(Projections[[f]])[1]
  mapStack <- stack()
  for(y in 1:numYears) {
    df <- data.frame(t(Projections[[f]][y,,,'Niche'] * sumWeight))
    df$ID <- row.names(df)
    map <- refRast
    map[as.integer(df$ID)] <- df$J+df$A
    names(map) <- paste0('year', y)
    mapStack <- addLayer(mapStack, map)
  }
  mapList[[f]] <- mapStack
}

numYears <- nlayers(mapList[[1]])
averages <- stack()
for(y in 1:numYears) {
  allMaps <- stack()
  for(r in projFiles) {
    allMaps <- addLayer(allMaps, mapList[[r]][[y]])
  }
  averageYear <- mean(allMaps)
  names(averageYear) <- paste0('year', y)
  averages <- addLayer(averages, averageYear)
}

writeRaster(averages, 'average_maps.img', format='HFA', overwrite=TRUE)

#### EOF ####
