#' @title Read jointly all output files from the network data event simulator.
#'
#' @description This function reads all output files from the network data event 
#' simulator from an output path name specified by the user and returns an rds 
#' file in another input path name for each of the files.
#'
#' @param path_output path name where the output csv files of the simulator are.
#' 
#' @return A named \code{list} with components:
#' 
#' \itemize{
#'  
#'  \item outputParameters
#'  \item antenna.dt
#'  \item events.dt
#'  \item gridParam
#'  \item RSS_Sim.dt
#'  \item coverArea
#'  \item coverArea_sf
#'  \item positionsConnections.dt
#'  \item anteAreas.dt.
#'  
#' }
#'
#' @examples
#' \dontrun{
#' path_output <- path_name1
#' get_simScenario(path_output)
#' 
#' }
#'
#' @import data.table
#'
#' @include tileEquivalence.R
#'
#' @export
get_simScenario <- function(fileGridName = NULL, 
                            filePersonsName = NULL,
                            fileEventsInfoName = NULL, 
                            fileAntennasPosName = NULL,
                            fileSignalName = NULL,
                            fileCoverName = NULL, 
                            antennasConfig.dt = NULL,
                            connectionType = NULL,
                            map){
  

  ### TechDebt: Under the assumption of only ONE MNO                        ####
  
  output <- list()
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  gridParFile                               #####
  
  if(!is.null(fileGridName)){
    
    cat('Reading and parsing gridParFile ...')
    gridParam <- fread(fileGridName, sep = ',', header = TRUE, stringsAsFactors = FALSE)
    ncol_grid  <- gridParam[['No Tiles Y']]
    nrow_grid  <- gridParam[['No Tiles X']]
    # tile_sizeX <- gridParam[['X Tile Dim']]
    # tile_sizeY <- gridParam[['Y Tile Dim']]
    # ntiles_x   <- gridParam[['No Tiles X']]
    # ntiles_y   <- gridParam[['No Tiles Y']]
    # ntiles     <- ntiles_x * ntiles_y
    
    cat('ok.\n')
    
    ### * Correspondence Grid tiles - R raster cells                          ####
    cat('     tile-raster correspondence...')
    tileEquiv.dt <- data.table(tileEquivalence(ncol_grid, nrow_grid))
    cat(' ok.\n')
    
    output$gridParam <- gridParam
    output$tileEquiv.dt <- tileEquiv.dt
  }
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  personsFile                               #####
  
  if(!is.null(filePersonsName)){
    
    cat('Reading and parsing persons.csv ...')
    persons.dt <- fread(filePersonsName, sep = '\n', stringsAsFactors = FALSE)
    colNames <- c('time', 'personID', 'x', 'y', 'tile', 'device_1', 'device_2')
    setnames(persons.dt, 'V1')
    persons_parsed.dt <- persons.dt[, tstrsplit(V1, split = ',')]
    setnames(persons_parsed.dt, colNames)
    persons_parsed.dt[
      , time := as.integer(time)][
        , x := as.numeric(x)][
          , y := as.numeric(y)][
            , tile := as.integer(tile)]
    #persons_parsed.dt <- merge(persons_parsed.dt, tileEquiv.dt, by = 'tile', all.x = TRUE)
    
    cat('ok.\n')

    output$persons_parsed.dt <- persons_parsed.dt
  }
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  eventsFile                                #####
  
  if(!is.null(fileEventsInfoName)){
    cat('Reading and parsing events.csv ...')
    
    allEvents.dt <- fread(fileEventsInfoName, sep = '\n', stringsAsFactors = FALSE, 
                          colClasses = c('character'))
    cols <- strsplit(names(allEvents.dt), split = ',')[[1]]
    allEvents.dt <- allEvents.dt[, tstrsplit(get(names(allEvents.dt)), split = ',')]
    setnames(allEvents.dt, cols)
    allEvents.dt <- allEvents.dt[!duplicated(allEvents.dt)]
    setnames(allEvents.dt , c('time', 'antennaID', 'eventCode', 'device', 'networkType', 'TA', 'x', 'y', 'tile'))
    allEvents.dt <- allEvents.dt[
      , antennaID := str_pad(antennaID, max(nchar(antennaID)), pad="0")][
      , time := as.integer(time)][
      , x := as.numeric(x)][
      , y := as.numeric(y)]
    
    allEvents.dt[, obsVar := do.call(paste, c(.SD, sep = "-")), .SDcols = c('antennaID', 'eventCode')]
    x <- allEvents.dt[, list(obsVar = paste(obsVar, collapse = ';')), by = c('device', 'time')]
    
    events.dt <- allEvents.dt[
      eventCode %in% c('0', '2', '3')]
    
    events.dt_noDup <- copy(events.dt)[, list(eventCode = as.character(min(as.numeric(eventCode)))), by = c("time", "device")]
    events.dt <- merge(events.dt_noDup, events.dt, by = names(events.dt_noDup), all.x = TRUE)
    events.dt <- events.dt[
      !duplicated(events.dt, by = c("time", "device", "eventCode"))][
      , .(time, device, eventCode, antennaID, obsVar)][
      order(time)]
    
    cat('ok.\n')
    
    output$events.dt <- events.dt

  }
    
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  antennaPosFile                            #####
  
  if(!is.null(fileAntennasPosName)){
    
    if(is.null(antennasConfig.dt)){
      stop("Error: antennasConfig.dt is missing.")
    }
    cat('Reading and parsing antennaPos.csv ...')
    
    antenna.dt<- fread(fileAntennasPosName, sep = ',')
    setnames(antenna.dt, c('time', 'antennaID', 'x', 'y', 'mnoID', 'tile'))
    antenna.dt <- antenna.dt[, x := NULL][, y := NULL]
    antenna.dt <- antenna.dt[, antennaID := str_pad(antennaID, max(nchar(antennaID)), pad="0")]
    antenna.dt <- antenna.dt[antennasConfig.dt, on = "antennaID"][, time := NULL]
    #antenna.dt <- merge(antenna.dt, tileEquiv.dt, by = 'tile', all.x = TRUE)
    cat('ok.\n')
    
    output$antenna.dt <- antenna.dt
  }
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  RSSfile                                   #####
  
  if(!is.null(fileSignalName)){
    
    if(is.null(connectionType)){
      stop("Error: connectionType is missing.")
    }
    if(is.null(antennasConfig.dt)){
      stop("Error: antennasConfig.dt is missing.")
    }
    cat('Reading and parsing RSS.csv ...')
    
    RSS_Sim.dt <- fread(fileSignalName, sep = ",", header = TRUE, 
                        stringsAsFactors = FALSE)
    nTiles <- dim(RSS_Sim.dt)[2] - 1
    setnames(RSS_Sim.dt, c('antennaID', 0:(nTiles - 1)))
    RSS_Sim.dt <- melt(RSS_Sim.dt, id.vars = 'antennaID', variable.name = 'tile', variable.factor = FALSE, value.name = 'RSS')
    RSS_Sim.dt <- RSS_Sim.dt[, antennaID := stringr::str_pad(antennaID, max(nchar(antennaID)), pad = "0")]
    # Transform signal strength to signal quality (SDM, signal dominance measure):
    RSS_Sim.dt <- merge(RSS_Sim.dt, antennasConfig.dt, all.x = TRUE, by = "antennaID")
    if(connectionType == "strength"){
      RSS_Sim.dt[, SDM := (1/(1 + exp(-SSteep * (RSS - Smid))))][
        , SDM_ori := SDM][
        , RSS_ori := RSS]
    }
    if(connectionType == "quality"){
      RSS_Sim.dt[, RSS_actual := (1/(SSteep * log(RSS/(1-RSS)) ))][
        , SDM_ori := RSS][
        , SDM := RSS][
        , RSS_ori := RSS_actual][
        , RSS := RSS_actual][
        , RSS_actual := NULL]
    }
    
    
    # Make RSS or SDM = NA if the tile is out the coverage area
    RSS_Sim.dt[
      , RSS := ifelse(RSS < Smin, NA, RSS)]
    if(any(names(RSS_Sim.dt) %in% "Qmin")){
      RSS_Sim.dt[
        , SDM := ifelse(SDM < Qmin, NA, SDM)]
    }
  
    RSS_Sim.dt[
      , tile := as.integer(tile)]
    
    minTile <- RSS_Sim.dt[, min(tile)]
    maxTile <- RSS_Sim.dt[, max(tile)]
      
    RSS_noNA.dt <- RSS_Sim.dt[!is.na(RSS)]
    nCover_tile.dt <- merge(RSS_noNA.dt[, .N, by = tile], data.table(tile = minTile:maxTile), all = TRUE)[
      is.na(N), N := 0]
    setnames(nCover_tile.dt, 'N', 'nAntCover')
    RSS_Sim.dt <- RSS_Sim.dt[nCover_tile.dt, on = 'tile']
    #RSS_Sim.dt <- merge(RSS_Sim.dt, tileEquiv.dt, by = 'tile', all.x = TRUE)
    
    removeCols <- c("power", "attenuationfactor", "Smin", "Qmin", "Smid",
                    "SSteep", "S0", "cover_radio")
    set(RSS_Sim.dt, j = removeCols, value = NULL)
  
    cat('ok.\n')
    
    output$RSS_Sim.dt <- RSS_Sim.dt
  }  
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      READ  coverageAreaFile                        #####
  
  if(!is.null(fileCoverName)){
    if(is.null(map)){
      stop("Error: map is missing.")
    }
    cat('Reading and parsing coverageArea.csv ...')

    coverArea <- fread(fileCoverName, sep ='\n', header = TRUE, stringsAsFactors = FALSE)
    setnames(coverArea, 'lines')
    coverArea[, antennaID := tstrsplit(lines, split = ',POLYGON')[[1]]]
    coverArea[, wkt := substring(lines, regexpr("POLYGON", lines))]
    coverArea <- coverArea[, c('antennaID', 'wkt'), with = FALSE]
    antennas <- coverArea[['antennaID']]
    antennas <- stringr::str_pad(antennas, max(nchar(antennas)), pad = "0")
    coverArea <- lapply(coverArea[['wkt']], function(wkt){
      
      polygon <- rgeos::readWKT(wkt)
      polygon <- raster::intersect(polygon, map)
      return(polygon)
    })
    names(coverArea) <- antennas
    coverArea_sf <- lapply(coverArea, sf::st_as_sf)
  
    cat('ok.\n')
    
    ### * Obtain area of coverage                                             ####
    cat('     anteAreas.dt ...')
    anteAreas <- sapply(coverArea, rgeos::gArea)
    anteAreas.dt <- data.table(antennaID = names(anteAreas),
                               cover_area = anteAreas)
    cat('ok.\n')
    
    output$coverArea <- coverArea
    output$coverArea_sf <- coverArea_sf
    output$anteAreas.dt <- anteAreas.dt
    
  }
  
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                      CREATE new data                                 #####
  

  if(!is.null(filePersonsName) & !is.null(fileEventsInfoName)){
    cat('Create new data ...\n')
    ### * Generate positionsConnections.dt                                    ####
    cat('     positionsConnections.dt ...')
      
    tempDT <- melt(persons_parsed.dt[!is.na(device_1)], id.vars = c('time', 'personID', 'x', 'y', 'tile'),
                   measure.vars = c('device_1', 'device_2'),
                   value.name = 'device')[
                   !is.na(device)][
                   , variable := NULL][
                  , time := as.integer(time)]

    positionsConnections.dt <- merge(
      events.dt, tempDT[,  .(time, device, personID, x, y, tile)],
      by = c("time", "device"), all = TRUE)[
      , time := as.integer(time)][
      , tile := as.integer(tile)]
    #positionsConnections.dt <- merge(
    #  positionsConnections.dt, tileEquiv.dt, by = 'tile', all.x = TRUE)
    rm(tempDT)
    cat('ok.\n')
    
    output$positionsConnections.dt <- positionsConnections.dt
  }
  

   
  ####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
  #####                           OUTPUT                                     #####
  

  return(output)
  
}
