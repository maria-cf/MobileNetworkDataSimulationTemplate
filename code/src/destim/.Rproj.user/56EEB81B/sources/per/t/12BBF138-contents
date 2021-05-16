####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
#path_root         <- 'C:/Users/U0xxxxx/Documents/MobileNetworkDataSimulationTemplate'
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

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)           # manage data
library(tmap)
library(mobvis)
library(mobloc)
library(xml2)
library(magrittr)
library(tibble)
library(tidyr)
library(stringr)
library(rgeos)
library(deduplication)
library(ggplot2)
library(latex2exp)

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
# Function to compute center of probabilities (to be merged with that from deduplication)
source(file.path(path_source, 'cp.R'))
# Function to compute radius of dispersion
source(file.path(path_source, 'rd.R'))


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                      READ RDS DATA FILES                               ####
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
tile_size  <- mean(c(tile_sizeX, tile_sizeY))
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

# Network event data
events.dt    <- simScenario.list$events.dt

# Centroid coordinates of each tile
centroidCoord.dt <- fread(file.path(path_resources, 'centroidCoord.csv'))

# True positions of each device
truePositions_device.dt <- fread(file.path(path_grTruth, 'truePositions_device.csv'))


# Posterior location probabilities
postLocProb.dt <- fread(file.path(path_postLoc, paste0('postLocProb_', geolocation_model, '_RSS-', geolocation_prior, '.csv')))
postLocProb.dt <- merge(postLocProb.dt, centroidCoord.dt, by = 'tile')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####      CENTER OF PROBABILITY AND ROOT MEAN SQUARED DISPERSION            ####

#### ** Bias                                                                ####
cp_truePos.dt <- postLocProb.dt[
  , as.list(cp(x = centroidCoord_x, y = centroidCoord_y, w = postLocProb)), by = c('device', 'time')][
  truePositions_device.dt[, .(time, device, x, y)], on = c('device', 'time')][
  , cp_tp := sqrt( (cp_x - x)**2 + (cp_y - y)**2)][
  , .(device, time, cp_tp)][
  order(device, time)]

#### ** rmsd                                                                ####
rmsd.dt <- postLocProb.dt[
  , list(rd = rd(x = centroidCoord_x, y = centroidCoord_y, w = postLocProb)), by = c('device', 'time')]

msd.dt <- merge(cp_truePos.dt, rmsd.dt, by = c('device', 'time'))
msd.dt[
  , msd := cp_tp**2 + rd**2][
  , model := paste0(geolocation_model, '-', geolocation_prior)]

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####              PLOT DISTRIBUTIONS of b, rmsd AND msd                     ####
#### ** Bias                                                                ####
ggplot(msd.dt, aes(x = model, y = cp_tp / tile_size)) +
  geom_boxplot(aes(fill = '')) +
  theme_bw() +
  labs(x = '', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(msd.dt, aes(x = time, y = cp_tp / tile_size, group = time)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(msd.dt, aes(x = device, y = cp_tp / tile_size, group = device)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nDevice', y = TeX('$b_{dt}$ (no. tiles)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 6, angle = 90),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


#### ** rmsd
ggplot(msd.dt, aes(x = model, y = rd / tile_size)) +
  geom_boxplot(aes(fill = model)) +
  theme_bw() +
  labs(x = '', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


ggplot(msd.dt, aes(x = time, y = rd / tile_size, group = time)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

ggplot(msd.dt, aes(x = device, y = rd / tile_size, group = device)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nDevice', y = 'rmsd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 6, angle = 90),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


#### ** msd
ggplot(msd.dt, aes(x = model, y = msd / tile_size**2)) +
  geom_boxplot(aes(fill = model)) +
  theme_bw() +
  labs(x = '', y = 'msd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


ggplot(msd.dt, aes(x = time, y = msd / tile_size**2, group = time)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nTime (s)', y = 'msd (no. tiles)\n', title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


ggplot(msd.dt, aes(x = device, y = msd / tile_size**2, group = device)) +
  geom_boxplot(aes(fill = model)) +
  facet_grid(model ~ . ) +
  theme_bw() +
  labs(x = '\nDevice', y = TeX('msd (no. tiles$^{2}$)\n'), title = '') +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 6, angle = 90),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')


####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####            PLOT DEVICE TIME SERIES OF of b, rmsd AND msd               ####
devID <- sort(unique(msd.dt$device))[5]
ggplot(msd.dt[device == devID], aes(x = time, y = cp_tp / tile_size, group = time)) +
  geom_point(aes(color = model)) +
  facet_grid(model ~ . ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = paste0('Device ', devID)) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

devID <- sort(unique(msd.dt$device))[5]
ggplot(msd.dt[device == devID], aes(x = time, y = rd / tile_size, group = time)) +
  geom_point(aes(color = model)) +
  facet_grid(model ~ . ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = paste0('Device ', devID)) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')

devID <- sort(unique(msd.dt$device))[5]
ggplot(msd.dt[device == devID], aes(x = time, y = msd / tile_size**2, group = time)) +
  geom_point(aes(color = model)) +
  facet_grid(model ~ . ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(x = '\nTime (s)', y = TeX('$b_{dt}$ (no. tiles)\n'), title = paste0('Device ', devID)) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.position = 'none')
