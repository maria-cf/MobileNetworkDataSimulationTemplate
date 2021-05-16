####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root         <- '/Users/Maria/Desktop/Root/MobileNetworkDataSimulationTemplate'
path_source       <- file.path(path_root, 'code/src')
path_simConfig    <- file.path(path_root, 'data/simulatorConfig')
path_events       <- file.path(path_root, 'data/networkEvents')
path_eventLoc     <- file.path(path_root, 'data/eventLocProb')
path_resources    <- file.path(path_root, 'param/resources')
path_processParam <- file.path(path_root, 'param/process')
path_postLoc      <- file.path(path_root, 'data/postLocProb')
path_img          <- file.path(path_root, 'metrics/img')
path_grTruth      <- file.path(path_root, 'param/groundTruth')
path_data         <- file.path(path_root, 'data')

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)           # manage data
library(tmap)
library(mobvis)
library(mobloc)
library(magrittr)
library(tibble)
library(latex2exp)
library(xml2)                 # to read simulatorConfig and parameters
library(tidyr)                # transform xml documents to tables (tibbles)
library(rgeos)                # spatial data
library(destim)               # HMM model
library(stringr)              # to pad strings
library(Matrix)

# Function get_simConfig to read the input files of the simulator
source(file.path(path_source, 'get_simConfig.R'))
# Function get_simScenario to read the output files of the simulator
source(file.path(path_source, 'get_simScenario.R'))
# Function tileEquivalence to compute the equivalence between rastercell (R) and tiles (simulator)
source(file.path(path_source, 'tileEquivalence.R'))


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####           LOAD PARAMETERS AND NETWORK EVENT DATA                     #####
simConfigParam.list <- get_simConfig(path_simConfig)

# grid
fileGridName       <- file.path(path_resources, simConfigParam.list$filesNames$fileGridName)
fileEventsInfoName <- file.path(path_events, simConfigParam.list$filesNames$fileEventsInfoName)
simScenario.list   <- get_simScenario(fileGridName = fileGridName, fileEventsInfoName = fileEventsInfoName)
gridParam  <- simScenario.list$gridParam
tile_sizeX <- gridParam[['X Tile Dim']]
tile_sizeY <- gridParam[['Y Tile Dim']]

# Network event data
events.dt    <- simScenario.list$events.dt
TA.dt <- fread(file.path(path_events, 'AntennaInfo_MNO_MNO1.csv'),
               colClasses = c('integer', 'character', 'character', 'character', 'character', 'character', rep('numeric', 3)))[
                 , .(t, PhoneId, NetworkType, TA)]
TA.dt <- TA.dt[!duplicated(TA.dt, by = c('t', 'PhoneId'))]
setnames(TA.dt, c('time', 'device', 'networkType', 'TA'))
events.dt <- merge(events.dt, TA.dt, by = c('time', 'device'))[
  , event := do.call(paste, c(.SD, sep = '_')), .SDcols = c('antennaID', 'TA')]

# time
t_increment <- simConfigParam.list$simConfigParameters$t_increment

# Maximum velocity
vMax_ms <- as.numeric(parameters.xml$process_parameters$geolocation$params$vmax_ms[[1]])

deviceIDs <- sort(unique(events.dt$device))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                      SET TIME PADDING PARAMETERS                       ####
distMax  <- vMax_ms * t_increment
pad_coef <- as.integer(ceiling(distMax / max(tile_sizeX, tile_sizeY)))
pad_coef <- pad_coef + 1

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                INITIAL STATE DISTRIBUTION (PRIOR)                      ####
initialDistr_RSS_network.dt <- RSS.dt[
  , watt := 10**( (RSS_ori - 30) / 10 )][
    , total := sum(watt, na.rm = TRUE)][
      , list(num = sum(watt, na.rm = TRUE), total = unique(total)), by = 'rasterCell'][
        , prior_network := num / total][
          order(rasterCell)]
prior_network <- initialDistr_RSS_network.dt$prior_network

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                          EMISSION MODEL                                ####
emissionProb_rasterCell.dt <- fread(file.path(path_eventLoc, "eventLocProb.csv"),
                                    colClasses = c('character', 'numeric', 'numeric', 'character', 'character', 'numeric'),
                                    sep = ',')[
                                      tileEquiv.dt, on = 'tile'][
                                        , c('device', 'time', 'tile') := NULL][
                                          , TA := as.integer(TA)][
                                            order(event_cellID, TA, rasterCell)]
setcolorder(emissionProb_rasterCell.dt, c('rasterCell', 'event_cellID', 'TA', 'eventLocProb'))
event_cellIDs <- sort(unique(emissionProb_rasterCell.dt$event_cellID))
TAs <- 0:max(emissionProb_rasterCell.dt$TA)
events_ID.dt <- as.data.table(expand.grid(event_cellID = event_cellIDs, TA = TAs))[
  order(event_cellID, TA)][
    , events := do.call(paste, c(.SD, sep = '_'))]
events_ID <- events_ID.dt$events
emissionProb_rasterCell.matrix <- as.matrix(
  dcast(emissionProb_rasterCell.dt,
        rasterCell ~ event_cellID + TA, value.var = 'eventLocProb')[
          , rasterCell := NULL])
newEvents <- setdiff(events_ID, dimnames(emissionProb_rasterCell.matrix)[[2]])
emissionProb_newEvents.matrix <- matrix(
  0, nrow = nrow(emissionProb_rasterCell.matrix),
  ncol = length(newEvents))
colnames(emissionProb_newEvents.matrix) <- newEvents
emissionProb_rasterCell.matrix <- cbind(
  emissionProb_rasterCell.matrix, emissionProb_newEvents.matrix)[, events_ID]
emissionProb_rasterCell.matrix[is.na(emissionProb_rasterCell.matrix)] <- 0

model <- HMMrectangle(nrow_grid, ncol_grid)
emissions(model) <- emissionProb_rasterCell.matrix # eventLoc for each antenna

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####                          TRANSITION MODEL                            #####
model <- initparams(model)  # initialize transition prob
model <- minparams(model)   # parameter reduction according to restrictions

istates(model) <- prior_network

isDiagonal <- vector('logical', length = length(deviceIDs))
names(isDiagonal) <- deviceIDs

for (devID in deviceIDs){

  events_device.dt <- events.dt[
    device == devID, .(device, time, event)][
      order(device, time)]

  events_deviceID  <- unlist(events_device.dt[, c("event")])

  events_deviceID_pad <- rep(NA, pad_coef * length(events_deviceID))
  events_deviceID_pad[seq(1, length(events_deviceID_pad), by = pad_coef)] <- events_deviceID

  colEvents <- sapply(events_deviceID_pad,
                      function(x) ifelse(!is.na(x),
                                         which (x == colnames(emissionProb_rasterCell.matrix)), NA))

  model_devID <- fit(model, colEvents, init = TRUE)

  mat <- model_devID$transitions
  mat[2, ] <- -1 * mat[2,]
  max_NonDiagonalTerm <- max(abs(model_devID$parameters$transitions[which(colSums(mat) != 0)]))
  isDiagonal[devID] <- (max_NonDiagonalTerm < 1e-7)

}

names(isDiagonal)[!isDiagonal] -> devicesID_nonDiagonal
diagonal <- data.frame(device=as.integer(deviceIDs),isDiagonal)

msd.dt <- readRDS(file.path(path_data,"msd.dt.rds"))
msd_nonDiagonal.dt <- msd.dt[as.character(device)%chin%devicesID_nonDiagonal]
msd.dt <- merge(msd.dt,diagonal,by = "device", all.x=TRUE)

#Con TA
msdTA.dt <- readRDS(file.path(path_data,"msd_TA.dt.rds"))
msdTA_nonDiagonal.dt <- msdTA.dt[as.character(device)%chin%devicesID_nonDiagonal]
msdTA.dt <- merge(msdTA.dt,diagonal,by = "device", all.x=TRUE)

