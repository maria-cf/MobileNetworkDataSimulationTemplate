####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root        <- 'C:/Users/U0xxxxx/Documents/MobileNetworkDataSimulationTemplate'
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

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                        SET PARAMETERS                                  ####
processParam.list <- as_list(read_xml(file.path(path_processParam, 'parameters_process.xml')))



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                COMPUTE EMISSION PROB MATRICES                          ####
RSS.dt <- fread(file.path(path_resources, 'RSS.csv'), 
                colClasses = c('numeric', 'character', 'numeric', 'character', 'integer',
                               'numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                               'integer', 'numeric'))
emissionModel <- 'RSS'
eventVarName <- emissionModel
emissionProb.dt <- compute_eventLocProb(RSS.dt, emissionModel, eventVarName, 
                                     by.grid = 'tile', by.event = 'antennaID')[
  order(antennaID, tile)]
setnames(emissionProb.dt, c('antennaID', 'eventLocProb'), c('event_cellID', 'eventLocProb'))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                    SAVE EMISSION MATRIX                                ####
fwrite(emissionProb.dt, file.path(path_eventLoc, "eventLocProb.csv"))
