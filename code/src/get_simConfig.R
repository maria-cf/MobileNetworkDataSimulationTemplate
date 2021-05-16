#' @title Read jointly all config files for the network data event simulator.
#'
#' @description This function reads all configuration files for the network data
#'  event simulator stored in the directory specified as input argument and 
#'  returns a list with the parsed parameters.
#'
#' @param path_simConfig path name where the config (xml and wkt) files are.
#' 
#' @param simulationFile_name filename of the input xml file for the simulation
#' components needed by the simulator.
#' 
#' @param antennasFile_name filename of the input xml file for the antenna 
#' components needed by the simulator.
#' 
#' @param mapFile_name filename of the input WKT file for the map component 
#' needed by the simulator.
#' 
#' @param personsFile_name filename of the input xml file for the persons 
#' components needed by the simulator.
#' 
#' @param probabilitiesFile_name filename of the input xml file for 
#' the probabilities component needed by the simulator.
#'
#' @return A named \code{list} with components simulation.xml, antennas.xml, 
#' map, persons.xml, and probabilities.xml with the contents of each input file,
#' respectively. The function automatically writes the contents of these files
#' (partially parsed) into a set of rds and RData files:
#' 
#' \itemize{
#'  
#'  \item filenames.RData with the names of all filenames.
#'  \item simConfigParameters.RData with the values of input parameters.
#'  \item simulation.xml.rds.
#'  \item antennasConfig.dt.rds
#'  \item map.rds
#'  
#' }
#'
#' @examples
#' \dontrun{
#' path_simConfig <- path_name
#' get_simConfig(path_simConfig)
#' 
#' }
#'
#' @import data.table rgeos
#'
#' @export
get_simConfig <- function(path_simConfig, 
                          simulationFile_name = "simulation.xml",
                          antennasFile_name = "antennas.xml",
                          mapFile_name = "map.wkt",
                          personsFile_name = "persons.xml",
                          probabilitiesFile_name = "probabilities.xml"){
  
 
  ### TechDebt: Under the assumption of only ONE MNO                        ####
  
   
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  simulationFile_name                       #####
  
  cat('Reading and parsing simulation.xml ...')
  simulation.xml <- xml2::as_list(xml2::read_xml(file.path(path_simConfig, simulationFile_name)))
  
  MNOName        <- simulation.xml$simulation$mno$mno_name[[1]]
  connectionType <- simulation.xml$simulation$connection_type[[1]]
  conn_threshold <- simulation.xml$simulation$conn_threshold[[1]]
  prob_secDev    <- as.numeric(simulation.xml$simulation$prob_sec_mobile_phone[[1]])
  
  initialTime <- as.integer(simulation.xml$simulation$start_time[[1]])
  finalTime   <- as.integer(simulation.xml$simulation$end_time[[1]])
  t_increment <- as.integer(simulation.xml$simulation$time_increment[[1]])
  times       <- seq(from = initialTime, to = (finalTime - t_increment), by = t_increment)
  
  fileGridName       <- simulation.xml$simulation$grid_file[[1]]
  filePersonsName    <- simulation.xml$simulation$persons_file[[1]]
  fileEventsInfoName <- paste0('AntennaInfo_MNO_', MNOName, '.csv')
  fileSignalName     <- paste0('SignalMeasure_', MNOName, '.csv')
  fileCoverName      <- paste0('AntennaCells_', MNOName, '.csv')
  
  cat('ok.\n')
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  antennasFile_name                         #####
  
  cat('Reading and parsing antennas.xml ...')
  antenna.xml <- as_list(read_xml(file.path(path_simConfig, antennasFile_name)))
  antenna.tib <- as_tibble(antenna.xml)%>%
    # new tidyr function
    unnest_wider(antennas) %>%
    # unnest same length list cols
    unnest(cols = names(.)) %>%
    unnest(cols = names(.)) %>%
    # convert using readr parser
    readr::type_convert()
  antennasConfig.dt <- as.data.table(antenna.tib)[
      , antennaID := as.character(.I)]
  
  keep <- intersect(names(antennasConfig.dt), 
                    c("antennaID", "power", "attenuationfactor", "Smid", 
                      "SSteep", "Smin", "Qmin", "x", "y", "mno_name", "maxconnections"))
  
  antennasConfig.dt <- antennasConfig.dt[, ..keep][
    , antennaID := str_pad(antennaID, max(nchar(antennaID)), pad="0")]
  antennasConfig.dt[, S0 := 30 + 10 * log10(power)][
    , cover_radio := 10^((S0 - Smin)/(10 * attenuationfactor))]
  cat('ok.\n')
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  mapFile_name                              #####
  
  cat('Reading and parsing map.wkt ...')
  map <- readWKT(readLines(file.path(path_simConfig, mapFile_name)))
  cat('ok.\n')
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  personsFile_name                    #####
  
  cat('Reading and parsing persons.xml ...')
  persons.xml  <- as_list(read_xml(file.path(path_simConfig, personsFile_name)))
  num_persons  <- as.numeric(persons.xml$persons$num_persons[[1]])
  min_age      <- as.numeric(persons.xml$persons$min_age[[1]])
  max_age      <- as.numeric(persons.xml$persons$max_age[[1]])
  speed_walk   <- as.numeric(persons.xml$persons$speed_walk[[1]])
  speed_car    <- as.numeric(persons.xml$persons$speed_car[[1]])
  percent_home <- as.numeric(persons.xml$persons$percent_home[[1]])
  cat('ok.\n')
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  probabilitiesFile_name                    #####
  
  cat('Reading and parsing probabilities.xml ...')
  probabilities.xml     <- as_list(read_xml(file.path(path_simConfig, probabilitiesFile_name)))
  probabilitiesName     <- probabilities.xml$probabilities$prob_file_name_prefix[[1]]
  fileProbabilitiesName <- paste0(probabilitiesName, '_', MNOName, '.csv')
  cat('ok.\n')
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                       OUTPUT                                         #####
  
  output <- list(
    simulation.xml = simulation.xml,
    antennas_xml = antennasConfig.dt[],
    map = map,
    persons.xml = persons.xml,
    probabilities.xml = probabilities.xml,
    filesNames = list(fileProbabilitiesName = fileProbabilitiesName, 
                      fileGridName = fileGridName, 
                      filePersonsName = filePersonsName, 
                      fileEventsInfoName = fileEventsInfoName, 
                      fileSignalName = fileSignalName, 
                      fileCoverName = fileCoverName),
    simConfigParameters = list(MNOName = MNOName, 
                           probabilitiesName = probabilitiesName, 
                           connectionType = connectionType, 
                           conn_threshold = conn_threshold, 
                           prob_secDev = prob_secDev,
                           initialTime = initialTime, 
                           finalTime = finalTime, 
                           t_increment = t_increment, 
                           times = times, 
                           num_persons = num_persons, 
                           min_age = min_age, 
                           max_age = max_age, 
                           speed_walk = speed_walk, 
                           speed_car = speed_car, 
                           percent_home = percent_home)
  )
  return(output)
}
