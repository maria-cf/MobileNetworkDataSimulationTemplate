####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root        <- 'C:/Users/U0xxxxx/Documents/MobileNetworkDataSimulationTemplate'
path_source      <- file.path(path_root, 'code/src')
path_simConfig   <- file.path(path_root, 'data/simulatorConfig')
path_events      <- file.path(path_root, 'data/networkEvents')
path_eventLoc    <- file.path(path_root, 'data/eventLocProb')
path_resources   <- file.path(path_root, 'param/resources')
path_groundTruth <- file.path(path_root, 'param/groundTruth')

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
####                     LOAD INPUT DATA AND PARAMETERS                     ####
# time
simConfigParam.list <- get_simConfig(path_simConfig)
initialTime <- as.integer(simConfigParam.list$simulation.xml$simulation$start_time)
finalTime   <- as.integer(simConfigParam.list$simulation.xml$simulation$end_time)
t_increment <- as.integer(simConfigParam.list$simulation.xml$simulation$time_increment)
times       <- seq(from = initialTime, to = (finalTime - t_increment), by = t_increment)


regions.dt <- fread(file.path(path_resources, 'regions.csv'))

centroidCoord.dt <- fread(file.path(path_resources, 'centroidCoord.csv'))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  GET SIMULATION CONFIG PARAMS                          ####
simConfigParam.list <- get_simConfig(path_simConfig)
connectionType <- simConfigParam.list$simConfigParameters$connectionType
antennasConfig.dt <- simConfigParam.list$antennas_xml



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  GET SIMULATION SCENARIO PARAMS                        ####
filePersonsName    <- file.path(path_groundTruth, simConfigParam.list$filesNames$filePersonsName)
fileEventsInfoName <- file.path(path_events, simConfigParam.list$filesNames$fileEventsInfoName)

simScenario.list <- get_simScenario(filePersonsName = filePersonsName,
                                    fileEventsInfoName = fileEventsInfoName,
                                    map = simConfigParam.list$map)


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                       BUILD TRUE POSITIONS                             ####
truePositions_person.dt <- simScenario.list$persons_parsed.dt
truePositions_person_region.dt <- merge(
  truePositions_person.dt, centroidCoord.dt, by = 'tile')

truePositions_device.dt <- simScenario.list$positionsConnections.dt
truePositions_device_region.dt <- merge(
  truePositions_device.dt, centroidCoord.dt, by = 'tile')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####     COMPUTE TRUE NUMBER OF INDIVIDUALS DETECTED BY MNO IN REGIONs      ####
Nnet_t_true_region.dt <- truePositions_device_region.dt[
  !duplicated(truePositions_device_region.dt, by = c('time', 'personID'))][
    , .N, by = c('time', 'region')]
Nnet_t_true_region.dt <- merge(data.table(expand.grid(time = times, 
                                                  region = sort(unique(regions.dt$region)))),
                               Nnet_t_true_region.dt, all = TRUE)[
  is.na(N), N := 0]
setnames(Nnet_t_true_region.dt, 'N', 'Nnet_true')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####        COMPUTE TRUE OD MATRIX OF INDIVIDUALS DETECTED BY MNO           ####
devIDs_Nnet <- sort(unique(truePositions_person_region.dt[!is.na(device_1)]$device_1))
od_to.dt <- truePositions_device_region.dt[
  device %in% devIDs_Nnet][
  , .(tile, time, device, region)][
  , time := time - t_increment][
  , time_to := time + t_increment][
  time >= 0]
setnames(od_to.dt, c('tile', 'region'), c('tile_to', 'region_to'))
od_from.dt <- truePositions_device_region.dt[
  device %in% devIDs_Nnet][
  , .(tile, time, device, region)][
  , time_from := time]
setnames(od_from.dt, c('tile', 'region'), c('tile_from', 'region_from'))
od_region.dt <- merge(od_from.dt, od_to.dt, by = c('device', 'time'))[
  , time := NULL][
  , .N, by = c('time_from', 'time_to', 'region_from', 'region_to')]
setnames(od_region.dt, 'N', 'Nnet_true')
NnetOD_true_region.dt <- merge(
  data.table(expand.grid(time_from   = times[-length(times)],
                         time_to     = times[-1],
                         region_from = unique(regions.dt$region), 
                         region_to   = unique(regions.dt$region))),
  od_region.dt, all.y = TRUE)[
    is.na(Nnet_true), Nnet_true := 0]
setcolorder(NnetOD_true_region.dt, c('time_from', 'time_to', 'region_from', 'region_to', 'Nnet_true'))


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####       COMPUTE TRUE NUMBER OF INDIVIDUALS IN TARGET POPULATION          ####
N_t_true_region.dt <- truePositions_person_region.dt[
  , .N, by = c('time', 'region')][
  order(time, region)]
setnames(N_t_true_region.dt, 'N', 'Ntarget_true')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####       COMPUTE TRUE OD MATRIX OF INDIVIDUALS IN TARGET POPULATION       ####
od_to.dt <- truePositions_person_region.dt[
  , .(tile, time, personID, region)][
  , time := time - t_increment][
  , time_to := time + t_increment][
  time >= 0]
setnames(od_to.dt, c('tile', 'region'), c('tile_to', 'region_to'))
od_from.dt <- truePositions_person_region.dt[
  , .(tile, time, personID, region)][
  , time_from := time]
setnames(od_from.dt, c('tile', 'region'), c('tile_from', 'region_from'))
od_region.dt <- merge(od_from.dt, od_to.dt, by = c('personID', 'time'))[
  , time := NULL][
  , .N, by = c('time_from', 'time_to', 'region_from', 'region_to')]
setnames(od_region.dt, 'N', 'Ntarget_true')
NOD_t_true_region.dt <- merge(
  data.table(expand.grid(time_from   = times[-length(times)],
                         time_to     = times[-1],
                         region_from = unique(regions.dt$region), 
                         region_to   = unique(regions.dt$region))),
  od_region.dt, all = TRUE)[
    is.na(Ntarget_true), Ntarget_true := 0]
setcolorder(NOD_t_true_region.dt, c('time_from', 'time_to', 'region_from', 'region_to', 'Ntarget_true'))



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                    COMPUTE TRUE NUMBER OF DEVICES                      ####
Ndev_t_true_region.dt <- truePositions_device_region.dt[
  , .N, by = c('time', 'region')][
    order(time, region)]
setnames(Ndev_t_true_region.dt, 'N', 'Ndev_true')

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                 COMPUTE TRUE REGION OD MATRIX OF DEVICES               ####
od_to.dt <- truePositions_device_region.dt[
  , .(tile, time, device, region)][
    , time := time - t_increment][
      , time_to := time + t_increment][
        time >= 0]
setnames(od_to.dt, c('tile', 'region'), c('tile_to', 'region_to'))
od_from.dt <- truePositions_device_region.dt[
  , .(tile, time, device, region)][
    , time_from := time]
setnames(od_from.dt, c('tile', 'region'), c('tile_from', 'region_from'))
od_region.dt <- merge(od_from.dt, od_to.dt, by = c('device', 'time'))[
  , time := NULL][
    , .N, by = c('time_from', 'time_to', 'region_from', 'region_to')]
setnames(od_region.dt, 'N', 'Ndev_true')

NdevOD_t_true_region.dt <- merge(
  data.table(expand.grid(time_from   = times[-length(times)],
                         time_to     = times[-1],
                         region_from = unique(regions.dt$region), 
                         region_to   = unique(regions.dt$region))),
  od_region.dt, all.y = TRUE)[
    is.na(Ndev_true), Ndev_true := 0]
setcolorder(NdevOD_t_true_region.dt, c('time_from', 'time_to', 'region_from', 'region_to', 'Ndev_true'))



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                        SAVE RESOURCES                                  ####
fwrite(truePositions_person.dt, file.path(path_groundTruth, "truePositions_person.csv"))
fwrite(truePositions_device.dt, file.path(path_groundTruth, "truePositions_device.csv"))

fwrite(Nnet_t_true_region.dt, file.path(path_groundTruth, "Nnet_true_region.csv"))
fwrite(NnetOD_true_region.dt, file.path(path_groundTruth, "NnetOD_true_region.csv"))

fwrite(Ndev_t_true_region.dt, file.path(path_groundTruth, "Ndev_true_region.csv"))
fwrite(NdevOD_t_true_region.dt, file.path(path_groundTruth, "NdevOD_true_region.csv"))

fwrite(N_t_true_region.dt, file.path(path_groundTruth, "Ntarget_true_region.csv"))
fwrite(NOD_t_true_region.dt, file.path(path_groundTruth, "NtargetOD_true_region.csv"))

