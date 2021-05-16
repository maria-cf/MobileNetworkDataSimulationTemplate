#' @title Fit a static posterior location probability model.
#'
#' @description This function computes the posterior location probabilities 
#' according to static models with priors \code{uniform} or \code{network}.
#'  
#' @param RSS.dt data.table with columns: antennaID, tile, RSS, SDM, eventLocProb.
#' 
#' @param uniform logical depending if the approach with uniform prior must be obtained, 
#' by default TRUE.
#' 
#' @param network logical depending if the approach with network prior must be obtained, 
#' by default TRUE.
#' 
#' @param eventLocModel character with options: strength or quality, depending on
#' the type of connection done.
#' 
#' @return A named \code{list} with components postLocProb_device.dt and 
#' postLocProb_device.dt with the key-value pair data.tables with these probabilities.
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
compute_staticModel <- function(
  RSS_Sim.dt, events.dt, prior, eventLocModel){  
  
  output <- list()
  
  if(prior == 'uniform'){
    #####                 Fit static model with uniform prior                 ####
    cat('       fitting static model with uniform prior...')
    postLocProb.dt <- copy(RSS_Sim.dt)[, c('antennaID', 'tile', 'eventLocProb'), with = FALSE]
    postLocProb.dt[, postLoc_uniform := eventLocProb / sum(eventLocProb), by = 'antennaID']
    postLocProb.dt[, eventLocProb := NULL]
    postLocProb.dt[, antennaID := as.character(antennaID)]
    tempDT <- dcast(postLocProb.dt, formula = antennaID ~ tile, value.var = 'postLoc_uniform')
    tempDT <- merge(tempDT, events.dt[, c('time', 'antennaID'), with = FALSE], by = 'antennaID')
    postLocProb_device.dt <- melt(tempDT, id.vars = c('antennaID', 'time'), 
                                            variable.name = 'tile', value.name = 'postLocProb', 
                                            variable.factor = FALSE, value.factor = FALSE)[
                                              , tile := as.integer(tile)][
                                              , device := devID][
                                              order(device, time, tile)]
    setcolorder(postLocProb_device.dt, c('device', 'time', 'tile', 'antennaID', 'postLocProb'))
     
    cat(' ok.\n')
    
    output$postLocProb <- postLocProb_device.dt
    
  }
  
  if(prior == 'network'){
  
    #####                 Fit static model with network prior                 ####
    cat('       fitting static model with network prior...')
    
    if(eventLocModel == "RSS"){
      postLocProb.dt <- copy(RSS_Sim.dt)[
        , num := sum(RSS, na.rm = TRUE), by = 'tile'][
          , prior_network := num / sum(RSS, na.rm = TRUE)]
    }
    if(eventLocModel == "SDM"){
      postLocProb.dt <- copy(RSS_Sim.dt)[
        , num := sum(SDM, na.rm = TRUE), by = 'tile'][
        , prior_network := num / sum(SDM, na.rm = TRUE)]
    }  
    postLocProb.dt <- postLocProb.dt[
      , unnormalPostLoc := eventLocProb * prior_network][
      , postLoc_network := unnormalPostLoc / sum(unnormalPostLoc), by = c('antennaID')]
    postLocProb.dt[, eventLocProb := NULL]
    postLocProb.dt[, antennaID := as.character(antennaID)]
    tempDT <- dcast(postLocProb.dt, formula = antennaID ~ tile, value.var = 'postLoc_network')
    tempDT <- merge(tempDT, events.dt[, c('time', 'antennaID'), with = FALSE], by = 'antennaID')
    postLocProb_device.dt <- melt(tempDT, id.vars = c('antennaID', 'time'), 
                                            variable.name = 'tile', value.name = 'postLocProb', 
                                            variable.factor = FALSE, value.factor = FALSE)[
                                              order(time)][
                                              , tile := as.integer(tile)][
                                              , device := devID][
                                              order(device, time, tile)]
    setcolorder(postLocProb_device.dt, c('device', 'time', 'tile', 'antennaID', 'postLocProb'))
    
    cat(' ok.\n')
    
    output$postLocProb <- postLocProb_device.dt
  }
  
  
  return(output)

}