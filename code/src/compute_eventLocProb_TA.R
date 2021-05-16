compute_eventLocProb_TA <- function(networkConfig.dt, TA.dt, model, eventVarName, by.grid = 'tile', by.event = 'antennaID'){

    tempDT.dt <- merge(networkConfig.dt, TA.dt, by = c(by.grid, by.event))[
      , watt := 10**( (get(eventVarName) - 30) / 10 )][
      , watt_norm := watt / sum(watt, na.rm = TRUE), by = by.grid][
      , area_norm := area / sum(area, na.rm = TRUE), by = c(by.grid, by.event)][
      , unnormProb := watt_norm * area_norm][  
      , eventLocProb := unnormProb / sum(unnormProb, na.rm = TRUE), by = by.grid][
      is.na(eventLocProb), eventLocProb := 0][
       , watt := NULL][
       , watt_norm := NULL][
       , area_norm := NULL][  
       , unnormProb := NULL]

  output.dt <- tempDT.dt[ 
    , c(by.event, 'TA', by.grid, 'eventLocProb'), with = FALSE][
    , device := NA_character_][
    , time := NA_real_][  
    order(get(by.event), TA, get(by.grid))]
  setcolorder(output.dt, c('device', 'time', by.grid, by.event, 'TA', 'eventLocProb'))
  
  return(output.dt)
  
}
