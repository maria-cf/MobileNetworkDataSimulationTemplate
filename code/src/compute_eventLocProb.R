 compute_eventLocProb <- function(networkConfig.dt, model, eventVarName, by.grid = 'tile', by.event = 'antennaID', as.matrix = FALSE){
  
  if (model == 'RSS') {
    
    tempDT.dt <- networkConfig.dt[
      , watt := 10**( (get(eventVarName) - 30) / 10 )][
      , eventLocProb := watt / sum(watt, na.rm = TRUE), by = by.grid][
      is.na(eventLocProb), eventLocProb := 0][
      , watt := NULL]
  
  }
  
  if (model == 'SDM') {
    
    tempDT.dt <- networkConfig.dt[
      , eventLocProb := get(eventVarName) / sum(get(eventVarName), na.rm = TRUE), by = by.grid][
        is.na(eventLocProb), eventLocProb := 0]
    
  }

  output.dt <- tempDT.dt[ 
    , c(by.event, by.grid, 'eventLocProb'), with = FALSE][
    , device := NA_character_][
    , time := NA_real_][  
    order(get(by.event), get(by.grid))]
  setcolorder(output.dt, c('device', 'time', by.grid, by.event, 'eventLocProb'))
    
  if (as.matrix) {
  
   frmla <- paste0(c(by.grid,  by.event), collapse = ' ~ ')  
   output.mat <- as.matrix(dcast(output.dt,  as.formula(frmla), value.var = 'eventLocProb')[, (by.grid) := NULL])
   dimnames(output.mat) <- list(
     as.character(unique(output.dt[[by.grid]])), 
     as.character(unique(output.dt[[by.event]])))
   return(output.mat)   
  }
  
  return(output.dt)
  
}
