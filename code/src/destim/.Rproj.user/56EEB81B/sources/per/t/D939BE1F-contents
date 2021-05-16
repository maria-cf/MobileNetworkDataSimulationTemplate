####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root        <- '/Users/Maria/Desktop/Root/MobileNetworkDataSimulationTemplate'
path_source       <- file.path(path_root, 'code/src')
path_simConfig    <- file.path(path_root, 'data/simulatorConfig')
path_events       <- file.path(path_root, 'data/networkEvents')
path_eventLoc     <- file.path(path_root, 'data/eventLocProb')
path_resources    <- file.path(path_root, 'param/resources')
path_processParam <- file.path(path_root, 'param/process')

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)
library(stringr)
library(xml2)       # to read simulatorConfig
library(tidyr)      # transform xml documents to tables (tibbles)
library(rgeos)      # spatial data


# Function computeEventLoc to create the emissions matrix
source(file.path(path_source, 'compute_eventLocProb.R'))
source(file.path(path_source, 'get_simConfig.R'))
source(file.path(path_source, 'compute_TAarea.R'))
source(file.path(path_source, 'compute_eventLocProb_TA.R'))


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                        SET PARAMETERS                                  ####
processParam.list   <- as_list(read_xml(file.path(path_processParam, 'parameters_process.xml')))
simConfigParam.list <- get_simConfig(path_simConfig)
fileGridName        <- file.path(path_resources, simConfigParam.list$filesNames$fileGridName)

connectionType      <- simConfigParam.list$simConfigParameters$connectionType
centroidCoord.dt    <- fread(file.path(path_resources, 'centroidCoord.csv'),
                             colClasses = c('numeric', 'numeric', rep('integer', 6)))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                COMPUTE EMISSION PROB MATRICES                          ####
RSS.dt            <- fread(file.path(path_resources, 'RSS.csv'),
                           colClasses = c('numeric', 'character', 'numeric',
                                          'character', 'integer', 'numeric',
                                          'numeric', 'numeric', 'numeric',
                                          'numeric', 'integer', 'numeric'))
emissionModel     <- processParam.list$process_parameters$geolocation$emission_model[[1]]
eventVarName      <- 'RSS_ori'

antennasConfig.dt <- simConfigParam.list$antennas_xml
grid.dt <- fread(fileGridName)

TAarea.dt <- compute_TAarea(antennasConfig.dt, grid.dt, '3G')
TAarea.dt <- merge(TAarea.dt,
                   centroidCoord.dt[, .(centroidCoord_x, centroidCoord_y, tile)],
                   by = c('centroidCoord_x', 'centroidCoord_y'))[
            , .(tile, antennaID, TA, area)]

emissionProb.dt <- compute_eventLocProb_TA(
  RSS.dt, TAarea.dt, emissionModel, eventVarName,
  by.grid = 'tile', by.event = 'antennaID')[
  order(tile, antennaID, TA)]
setnames(emissionProb.dt, c('antennaID', 'eventLocProb'), c('event_cellID', 'eventLocProb'))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                    SAVE EMISSION MATRIX                                ####
fwrite(emissionProb.dt, file.path(path_eventLoc, "eventLocProb.csv"))
