####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root      <- 'C:/Users/U0xxxxx/Documents/MobileNetworkDataSimulationTemplate'
path_source    <- file.path(path_root, 'code/src')
path_simConfig <- file.path(path_root, 'data/simulatorConfig')
path_events    <- file.path(path_root, 'data/networkEvents')
path_eventLoc  <- file.path(path_root, 'data/eventLocProb')
path_resources <- file.path(path_root, 'param/resources')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)
library(stringr)
library(xml2)       # to read simulatorConfig
library(tidyr)      # transform xml documents to tables (tibbles)
library(rgeos)      # spatial data


# Function get_simConfig to read the input files of the simulator
source(file.path(path_source, 'get_simConfig.R'))
# Function get_simScenario to read the output files of the simulator
source(file.path(path_source, 'get_simScenario.R'))
# Function tileEquivalence to compute the equivalence between rastercell (R) and tiles (simulator)
source(file.path(path_source, 'tileEquivalence.R'))
# Function build_centroid_regions to create the data.table with the centroids and the regions
source(file.path(path_source, 'build_centroid_regions.R'))



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  GET SIMULATION CONFIG PARAMS                          ####
simConfigParam.list <- get_simConfig(path_simConfig)
connectionType      <- simConfigParam.list$simConfigParameters$connectionType
antennasConfig.dt   <- simConfigParam.list$antennas_xml



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  GET SIMULATION SCENARIO PARAMS                        ####
fileSignalName   <- file.path(path_resources, simConfigParam.list$filesNames$fileSignalName)
fileGridName     <- file.path(path_resources, simConfigParam.list$filesNames$fileGridName)

simScenario.list <- get_simScenario(fileSignalName = fileSignalName,
                                    fileGridName = fileGridName,
                                    antennasConfig.dt = antennasConfig.dt,
                                    connectionType = connectionType)

gridParam  <- simScenario.list$gridParam
ncol_grid  <- gridParam[['No Tiles Y']]
nrow_grid  <- gridParam[['No Tiles X']]
tile_sizeX <- gridParam[['X Tile Dim']]
tile_sizeY <- gridParam[['Y Tile Dim']]
ntiles_x   <- gridParam[['No Tiles X']]
ntiles_y   <- gridParam[['No Tiles Y']]


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  BUILD TILE/RASTER EQUIVALENCE                         ####
tileEquiv.dt <- simScenario.list$tileEquiv.dt


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                    BUILD SIGNAL DATA TABLE                             ####
RSS.dt <- simScenario.list$RSS_Sim.dt
RSS.dt <- merge(RSS.dt, tileEquiv.dt, by = 'tile', all.x = TRUE)


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                    BUILD CENTROID REGIONS                              ####
centroidCoord.dt <- build_centroid_regions(
  ntiles_x, ntiles_y, tile_sizeX, tile_sizeY, tileEquiv.dt)


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                        SAVE RESOURCES                                  ####
fwrite(RSS.dt, file.path(path_resources, "RSS.csv"))

fwrite(tileEquiv.dt, file.path(path_resources, 'tileEquiv.csv'))

fwrite(centroidCoord.dt, file.path(path_resources, "centroidCoord.csv"))

fwrite(centroidCoord.dt[, .(tile, region)], 
       file.path(path_resources, 'regions.csv'),
       row.names = FALSE, col.names = TRUE, sep = ',')

