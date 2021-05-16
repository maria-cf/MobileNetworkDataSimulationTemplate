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


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)           # manage data
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
# Function to fit and compute HMM model with the events of a specific device
source(file.path(path_source, 'compute_HMM.R'))
# Function to transform de output of compute_HMM
source(file.path(path_source, 'transform_postLoc.R'))
# Function to fit static model with uniform and network priors
source(file.path(path_source, 'compute_staticModel.R'))



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####           LOAD PARAMETERS AND NETWORK EVENT DATA                     #####
parameters.xml    <- as_list(read_xml(file.path(path_processParam, "parameters_process.xml")))
geolocation_model <- parameters.xml$process_parameters$geolocation$model_name[[1]]
geolocation_prior <- parameters.xml$process_parameters$geolocation$prior[[1]]
emission_model    <- 'RSS'
transition_model  <- parameters.xml$process_parameters$geolocation$transition_model[[1]]

simConfigParam.list <- get_simConfig(path_simConfig)

# time
times       <- simConfigParam.list$simConfigParameters$times
t_increment <- simConfigParam.list$simConfigParameters$t_increment

# grid
fileGridName       <- file.path(path_resources, simConfigParam.list$filesNames$fileGridName)
fileEventsInfoName <- file.path(path_events, simConfigParam.list$filesNames$fileEventsInfoName)
simScenario.list   <- get_simScenario(fileGridName = fileGridName, fileEventsInfoName = fileEventsInfoName)

gridParam  <- simScenario.list$gridParam
ncol_grid  <- gridParam[['No Tiles Y']]
nrow_grid  <- gridParam[['No Tiles X']]
tile_sizeX <- gridParam[['X Tile Dim']]
tile_sizeY <- gridParam[['Y Tile Dim']]
ntiles_x   <- gridParam[['No Tiles X']]
ntiles_y   <- gridParam[['No Tiles Y']]
ntiles     <- ntiles_x * ntiles_y

# tile-rasterCell equivalence
tileEquiv.dt <- simScenario.list$tileEquiv.dt

# RSS and other network parameters
RSS.dt <- fread(file.path(path_resources, 'RSS.csv'),
                colClasses = c('numeric', 'character', 'numeric', 'character', 'integer',
                               'numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                               'integer', 'numeric'))

# Maximum velocity
vMax_ms <- as.numeric(parameters.xml$process_parameters$geolocation$params$vmax_ms[[1]])

# Network event data
events.dt    <- simScenario.list$events.dt


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                      SET TIME PADDING PARAMETERS                       ####
distMax  <- vMax_ms * t_increment
pad_coef <- as.integer(ceiling(distMax / max(tile_sizeX, tile_sizeY)))
pad_coef <- pad_coef + 1



####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                INITIAL STATE DISTRIBUTION (PRIOR)                      ####
# Prepare prior_network distribution

if(geolocation_prior == "uniform"){

  prior_network <- rep(1 / ntiles, ntiles)

}

if(geolocation_prior == "network"){

  if(emission_model == "RSS"){

    initialDistr_RSS_network.dt <- RSS.dt[
      , watt := 10**( (RSS_ori - 30) / 10 )][
      , total := sum(watt, na.rm = TRUE)][
      , list(num = sum(watt, na.rm = TRUE), total = unique(total)), by = 'rasterCell'][
      , prior_network := num / total][
      order(rasterCell)]
    prior_network <- initialDistr_RSS_network.dt$prior_network

  }

  if(emission_model == "SDM"){

    initialDistr_SDM_network.dt <- RSS.dt[
      , total := sum(SDM_ori, na.rm = TRUE)][
      , list(num = sum(SDM_ori, na.rm = TRUE), total = unique(total)), by = 'rasterCell'][
      , prior_network := num / total][
      order(rasterCell)]
    prior_network <- initialDistr_SDM_network.dt$prior_network

  }
}

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                          EMISSION MODEL                                ####
# Load the event location probabilities:
# (i)  For the HMM this matrix contains the emission probabilities
# (ii) For the static analyses this matrix contains the likelihoods to apply Bayes' rule
emissionProb_rasterCell.dt <- fread(file.path(path_eventLoc, "eventLocProb.csv"),
                                    colClasses = c('character', 'numeric', 'numeric', 'character', 'numeric'),
                                    sep = ',')[
  tileEquiv.dt, on = 'tile'][
  , c('device', 'time', 'tile') := NULL][
  order(event_cellID, rasterCell)]
setcolorder(emissionProb_rasterCell.dt, c('rasterCell', 'event_cellID', 'eventLocProb'))
emissionProb_rasterCell.matrix <- as.matrix(
  dcast(emissionProb_rasterCell.dt,
        rasterCell ~ event_cellID, value.var = 'eventLocProb')[
  , rasterCell := NULL])


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####                          TRANSITION MODEL                            #####

if(transition_model == "HMMrectangle"){

  ## Initialize HMM                                                           ####
  cat('Initializing HMM...')
  model <- HMMrectangle(nrow_grid, ncol_grid)
  emissions(model) <- emissionProb_rasterCell.matrix # eventLoc for each antenna

  model <- initparams(model)  # initialize transition prob
  model <- minparams(model)   # parameter reduction according to restrictions

  istates(model) <- prior_network

  cat('ok.\n\n')

}


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####     COMPUTE POSTERIOR LOCATION PROBABILITIES AND INDICATORS           #####
cat('Computing posterior location probabilities...\n')

postLocProb.list      <- list()
postLocJointProb.list <- list()

deviceIDs <- sort(unique(events.dt$device))
RSS.dt <- merge(RSS.dt, emissionProb_rasterCell.dt[, .(rasterCell, event_cellID, eventLocProb)],
                by.x = c('rasterCell', 'antennaID'), by.y = c('rasterCell', 'event_cellID'),
                all.x = TRUE)
for(i in seq(along = deviceIDs)){

  ###                  :: For each device                          ####
  devID <- deviceIDs[i]
  cat(paste0('    device ', devID,'...\n'))

  # Selecting device events
  cat('       selecting network events...')
  events_device.dt <- events.dt[
    device == devID, .(device, time, antennaID)][
    order(device, time)]

  antennas_deviceID  <- unlist(events_device.dt[, c("antennaID")])
  cat(' ok.\n')

  if(!all(is.na(antennas_deviceID))){

    if(geolocation_model == "static"){

      if(geolocation_prior == "uniform"){

        #####                 Fit static model with uniform and network prior   ####
        cat('        fit static model with uniform prior...')
        staticModel_output <- compute_staticModel(RSS.dt, events_device.dt, prior = 'uniform')

        postLocProb_device.dt <- staticModel_output$postLocProb[
          tileEquiv.dt, on = 'tile']  ## ¿¿PARA QUÉ HACEMOS ESTO??
        postLocProb.list[[devID]] <- postLocProb_device.dt[
          postLocProb < 0, postLocProb := 0][
          , postLocProb := postLocProb / sum(postLocProb), by = 'time']
        fwrite(postLocProb_device.dt[, .(device, time, tile, antennaID, postLocProb)],
                    file.path(path_postLoc,
                              paste0('postLocProb_', geolocation_model, '_', emission_model,
                                     '_', geolocation_prior, '_', devID, '.csv')),
                    col.names = FALSE, row.names = FALSE, sep = ',')
      }

      if(geolocation_prior == "network"){

        #####                 Fit static model with network prior   ####
        cat('        fit static model with network prior...')
        staticModel_output <- compute_staticModel(RSS.dt, events_device.dt, prior = 'uniform')

        postLocProb_device.dt <- staticModel_output$postLocProb[
          tileEquiv.dt, on = 'tile']
        postLocProb.list[[devID]] <- postLocProb_device.dt[
          postLocProb < 0, postLocProb := 0][
          , postLocProb := postLocProb / sum(postLocProb), by = 'time']
        fwrite(postLocProb_device.dt[, .(device, time, tile, antennaID, postLocProb)],
                    file.path(path_postLoc,
                              paste0('postLocProb_', geolocation_model, '_', emission_model,
                                     '_', geolocation_prior, '_', devID, '.csv')),
                    col.names = FALSE, row.names = FALSE, sep = ',')

      }


    } # end if static


    if(geolocation_model == "HMM"){

      #####                Fit and compute HMM model                 ####
      cat('       fit and compute HMM model...\n')
      HMMmodel_output <- compute_HMM(model = model, observedValues = antennas_deviceID,
                                     pad_coef = pad_coef, init = TRUE)

      # model_devID <- HMMmodel_output$model_devID
      postLocProb_HMM_deviceID.matrix      <- HMMmodel_output$postLocP
      postLocJointProb_HMM_deviceID.matrix <- HMMmodel_output$postJointLocP

      cat(' ok.\n')
      #####                 Transform output HMM model                 ####
      cat('       transform output HMM model...\n')
      transform_output <- transform_postLoc(postLocP = postLocProb_HMM_deviceID.matrix,
                                           postLocJointP = postLocJointProb_HMM_deviceID.matrix,
                                           observedValues = antennas_deviceID,
                                           times = times, t_increment = t_increment,
                                           ntiles = ntiles, pad_coef = pad_coef,
                                           tileEquiv.dt = tileEquiv.dt, devID = devID,
                                           sparse_postLocP = TRUE, sparse_postLocJointP = TRUE)
      rm(postLocProb_HMM_deviceID.matrix)
      rm(postLocJointProb_HMM_deviceID.matrix)
      gc()
      postLocProb_deviceID.dt <- transform_output$postLocProb
      postLocProb.list[[as.character(devID)]] <- postLocProb_deviceID.dt[
        , c('device', 'time', 'tile', 'event_cellID', 'postLocProb'), with = FALSE]

      postLocJointProb_deviceID.dt <- transform_output$postLocJointProb
      postLocJointProb.list[[as.character(devID)]] <- postLocJointProb_deviceID.dt[
        , time_to := time_from + t_increment][
        , c('device', 'time_from', 'time_to', 'tile_from', 'tile_to', 'event_cellID_from', 'event_cellID_to', 'postLocProb'), with = FALSE]

      fwrite(postLocProb_deviceID.dt[, .(device, time, tile, event_cellID, postLocProb)],
             file.path(path_postLoc,
                       paste0('postLocProb_', geolocation_model, '_', emission_model,
                              '-', geolocation_prior, '_', devID, '.csv')),
             col.names = FALSE, row.names = FALSE, sep = ',')

      fwrite(postLocJointProb_deviceID.dt[, .(device, time_from, time_to, tile_from, tile_to, event_cellID_from, event_cellID_to, postLocProb)],
             file.path(path_postLoc,
                       paste0('postLocJointProb_', geolocation_model, '_', emission_model,
                              '-', geolocation_prior, '_', devID, '.csv')),
              col.names = FALSE, row.names = FALSE, sep = ',')

      rm(postLocJointProb_deviceID.dt)
      rm(transform_output)
      gc()
      cat(' ok.\n')
    }

  }
  ###                  :: End For each device                          ####

}


postLocProb.dt      <- rbindlist(postLocProb.list)
if(geolocation_model == "HMM"){
  postLocJointProb.dt <- rbindlist(postLocJointProb.list)
}

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####                   SAVE POSTERIOR LOCATION PROBS                      #####
fwrite(postLocProb.dt,
       file.path(path_postLoc,
                 paste0('postLocProb_', geolocation_model, '_',
                        emission_model, '-', geolocation_prior, '.csv')))

if(geolocation_model == "HMM"){

  fwrite(postLocJointProb.dt,
         file.path(path_postLoc,
                   paste0('postLocJointProb_', geolocation_model, '_',
                          emission_model, '-', geolocation_prior, '.csv')))
}
