#' @title Transform a location probability matrix into a key-value pair.
#'
#' @description This function transforms the matrix output for posterior 
#' location probabilities for each time instant (tiles $\times$ times) into a
#' key-value pair data.table with key = c(device, time, tile).
#' 
#' The function does also this transformation upon the matrix of joint posterior
#' location probabilities.
#' 
#' @param postLocP matrix with dimensions tiles $\times$ times with the 
#' posterior location probabilities.
#' 
#' @param postLocJointProb postLocP matrix of dimensions tiles $\times$ tiles * (times - 1) 
#' with the posterior joint location probabilities.
#' 
#' @param observedValues Vector with the observed values of the event variables.
#' 
#' @param times vector of time instants.
#' 
#' @param t_increment value with the increment of time between events.
#' 
#' @param ntiles total number of tiles.
#' 
#' @param pad_coef Integer coefficient denoting the number of time periods 
#' needed for the time padding.
#' 
#' @param tileEquiv.dt data.table with the equivalence between tiles and raster
#' cells according to the network event data simulator and rasters in R.
#' 
#' @para devID device ID.
#'
#' @return A named \code{list} with components postLocProb_HMM_deviceID.dt and 
#' postLocJointProb with the key-value pair data.tables with these probabilities.
#'  
#' }
#'
#' @examples
#' \dontrun{
#' path_input <- path_name1
#' path_rds   <- path_name2
#' read_inputSimulator(path_input, path_rds)
#' 
#' }
#'
#' @import data.table
#'
#' @include tileEquivalence.R
#'
#' @export
transform_postLoc <- function(postLocP, postLocJointP = NULL, observedValues,
                              times, t_increment, ntiles, pad_coef = 1, tileEquiv.dt, devID,
                              sparse_postLocP = FALSE, sparse_postLocJointP = FALSE){
  
  #### Checks over the input arguments                                                          ####
  
  # Check if postLocP
  # if () {
  #   
  #   stop(paste0('[transform_postLoc] .\n'))
  # }
  # 
  
  cat('Formatting output data.tables...\n')
  
  cat('         postLocProb...')
  postLocP0 <- copy(postLocP)
  postLocP0 <- as(postLocP0, "dgTMatrix")
  postLocP.dt <- data.table(rasterCell = postLocP0@i+1, 
                            time = postLocP0@j+1, 
                            postLocProb = postLocP0@x)
  # postLocP <- postLocP[, rasterCell := c(1:nrow(postLocP))]
  # postLocP <- melt(postLocP, id.vars = c('rasterCell'), variable.name = 'time', variable.factor = FALSE, value.name = 'postLocProb')[
  #   , time := as.integer(substr(time, 2, nchar(time)))]
  if(sparse_postLocP){
    
    postLocP.dt <- postLocP.dt[
      tileEquiv.dt, on = 'rasterCell'][
      , time := as.integer(time)]
  
  }

  if(!sparse_postLocP){
    
    allTiles.dt <- data.table(rasterCell = c(1:ntiles))
    allTiles.dt[, aux := 1]
    times.dt <- data.table(time = sort(unique(postLocP.dt$time)), aux = 1)
    allTilesTimes <- merge(allTiles.dt, times.dt, by = "aux", allow.cartesian = TRUE)
    allTilesTimes[, aux := NULL]
    allTilesTimes <- allTilesTimes[
      tileEquiv.dt, on = 'rasterCell'][
      , time := as.integer(time)]

    postLocP.dt <- merge(allTilesTimes, postLocP.dt, 
                         all.x = TRUE, by = intersect(names(allTilesTimes), names(postLocP.dt)))
    postLocP.dt[is.na(postLocProb), postLocProb := 0]
  
    }

  postLocP.dt <- postLocP.dt[time %in% seq(1, length(observedValues) * pad_coef, by = pad_coef)]
  names(times) <- sort(unique(postLocP.dt$time))
  postLocP.dt[, time := times[as.character(time)]][
    , device := devID]
  postLocP.dt <- postLocP.dt[
    data.table(time = times, event_cellID = observedValues), on = 'time']
  setcolorder(postLocP.dt, c('device', 'time', 'tile', 'rasterCell', 'event_cellID', 'postLocProb'))
  cat('      ok.\n')
  output <- list(postLocProb = postLocP.dt, postLocJointProb = NULL)
  rm(postLocP0)
  rm(postLocP.dt)

  if(!is.null(postLocJointP)){
    
    cat('         postLocJointProb...')
#    t1 <- Sys.time()
    postLocJointProb0 <- copy(postLocJointP)
    postLocJointProb0 <- as(postLocJointProb0, "dgTMatrix")
    postLocJointProb.dt <- data.table(rasterCell_from = postLocJointProb0@i+1, 
                                   j = postLocJointProb0@j+1, 
                                   postLocProb = postLocJointProb0@x)
    postLocJointProb.dt <- postLocJointProb.dt[, time_from := floor((j - 1) / ntiles) + 1][, rasterCell_to := j - ((time_from - 1) * ntiles)]
    postLocJointProb.dt <- postLocJointProb.dt[, j := NULL]
  
    postLocJointProb.dt <- postLocJointProb.dt[time_from %in% seq(1, length(observedValues) * pad_coef, by = pad_coef)]
    names(times) <- sort(unique(postLocJointProb.dt$time_from))
    postLocJointProb.dt[, time_from := times[as.character(time_from)]][
      , device := devID]
    postLocJointProb.dt <- merge(postLocJointProb.dt, tileEquiv.dt, by.x = 'rasterCell_from', by.y = 'rasterCell')
    setnames(postLocJointProb.dt, 'tile', 'tile_from')
    postLocJointProb.dt <- merge(postLocJointProb.dt, tileEquiv.dt, by.x = 'rasterCell_to', by.y = 'rasterCell')
    setnames(postLocJointProb.dt, 'tile', 'tile_to')

    if(!sparse_postLocJointP){
#      t2 <- Sys.time()
      allTiles.dt <- data.table(expand.grid(c(1:ntiles), c(1:ntiles)))
      setnames(allTiles.dt, c("Var1", "Var2"), c("rasterCell_from", "rasterCell_to"))
      allTiles.dt[, aux := 1]
      times.dt <- data.table(time_from = sort(unique(postLocJointProb.dt$time_from)), aux = 1)
      allTilesTimes <- merge(allTiles.dt, times.dt, by = "aux", allow.cartesian = TRUE)
      allTilesTimes[, aux := NULL][, device := devID]
      allTilesTimes <- merge(allTilesTimes, tileEquiv.dt, by.x = 'rasterCell_from', by.y = 'rasterCell')
      setnames(allTilesTimes, 'tile', 'tile_from')
      allTilesTimes <- merge(allTilesTimes, tileEquiv.dt, by.x = 'rasterCell_to', by.y = 'rasterCell')
      setnames(allTilesTimes, 'tile', 'tile_to')
      
      postLocJointProb.dt <- merge(allTilesTimes, postLocJointProb.dt, 
                                all.x = TRUE, by = intersect(names(allTilesTimes), names(postLocJointProb.dt)))
      rm(allTilesTimes)
      postLocJointProb.dt[is.na(postLocProb), postLocProb := 0]
      # t3 <- Sys.time()
      # difftime(t2, t1, units = "secs")
      # difftime(t3, t2, units = "secs")
    }

    setcolorder(postLocJointProb.dt, 
                c('device', 'time_from', 'tile_from', 'tile_to', 
                  'rasterCell_from', 'rasterCell_to', 'postLocProb'))
    postLocJointProb.dt <- postLocJointProb.dt[
      data.table(time_from = times[-length(times)]), on = 'time_from'][
      data.table(time_from = times[-length(times)], event_cellID_from = observedValues[-length(observedValues)]), on = 'time_from'][
      data.table(time_from = times[-length(times)], event_cellID_to = observedValues[-1]), on = 'time_from']

    output$postLocJointProb <- postLocJointProb.dt
    rm(postLocJointProb0)
    rm(postLocJointProb.dt)
    gc()
    cat(' ok.\n')
  }
  cat(' OK.\n')
  
  return(output)

}