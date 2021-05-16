####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
path_root         <- 'C:/Users/U0xxxxx/Documents/MobileNetworkDataSimulationTemplate'
path_source       <- file.path(path_root, 'code/src')
path_simConfig    <- file.path(path_root, 'data/simulatorConfig')
path_events       <- file.path(path_root, 'data/networkEvents')
path_eventLoc     <- file.path(path_root, 'data/eventLocProb')
path_resources    <- file.path(path_root, 'param/resources')
path_processParam <- file.path(path_root, 'param/process')
path_postLoc      <- file.path(path_root, 'data/postLocProb')
path_img          <- file.path(path_root, 'metrics/img')

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(data.table)           # manage data
library(tmap)
library(mobvis)
library(mobloc)

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
#####                 LOAD PARAMETERS AND NETWORK DAT                      #####
tileEquivalence.dt <- fread(file.path(path_resources, 'tileEquiv.csv'))
sim <- list(simConfig_dir = path_simConfig,
            resources_dir = path_resources,
            networkEvents_dir = path_events,
            postLocProb_dir = path_postLoc,
            mno = "MNO1",
            crs = sf::st_crs(NA))

# device to plot
devID <- '606'


rst      <- sim_get_raster(sim)
cp       <- sim_get_cellplan(sim)
map      <- sim_get_region(sim)
strength <- sim_get_signal_strength(sim, rst, cp)[
  , dist:= NA]
traj     <- sim_get_trajectory_data_TA(sim, device = devID)


param        <- mobloc_param()
strength_llh <- create_strength_llh(strength, param = param)


postLocProb <- sim_get_prob(sim, device = devID, 'postLocProb_HMM_RSS_TA-network')
setnames(postLocProb, c('tile'), c('rid'))
setcolorder(postLocProb, c('t', 'dev', 'rid', 'cell', 'p'))

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
#####           PLOT POSTERIOR LOCATION PROBABILITIES                      #####
# Retrieve and change the visualization settings
settings_plot <- mobvis_settings(cell_size = 2, dev_size = 2)
settings_plot$titles["pga"] <- "Posterior probabilities"


# Create the animation
settings_anim <- mobvis_settings_animation()
settings_anim$dev_color <- 'blue'
settings_anim$dev_size <- 5
animate_p(rst = rst,
          dt = postLocProb,
          cp = cp,
          traj = traj,
          region = map,
          settings = settings_anim,
          title = "Posterior distribution at time %s",
          filename = file.path(path_img, "postLocProb_dev_%s.mp4"),
          width = 700,
          height = 700,
          fps = 3)



