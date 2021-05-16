compute_TAarea <- function(antennaConf, grid, network_type){
  
  if (network_type == '4G') {
    
    delta <- 78.07 
    
  }
  
  if (network_type == '3G') {
    
    delta <- 554 
  }
  
  tileDim_x <- grid$'X Tile Dim'
  tileDim_y <- grid$'Y Tile Dim'
  nTiles_x  <- grid$'No Tiles X'
  nTiles_y  <- grid$'No Tiles Y'
  antennaConf <- as.data.table(antennaConf)
  nAnt <- nrow(antennaConf)
  output.dt <- data.table(antennaID = character(0), 
    centroidCoord_x = numeric(0), centroidCoord_y = numeric(0), 
    TA = integer(0), area = numeric(0))

  for (k in 1:nAnt) {
    cat(paste0('Antenna ', k, '...'))
    xa <- antennaConf[k]$x
    ya <- antennaConf[k]$y
    for (i in 1:nTiles_x) {
      for (j in 1:nTiles_y) {
        xc <- (0.5 + (i - 1)) * tileDim_x
        yc <- (0.5 + (j - 1)) * tileDim_y
        areaParcial <- 0
        areaTotal   <- 0
        left  <- xc - tileDim_x/2
        up    <- yc + tileDim_y/2
        right <- xc + tileDim_x/2
        down  <- yc - tileDim_y/2
        TAmin <- trunc(min(sqrt((xa-left)^2+(ya-up)^2),
                           sqrt((xa-right)^2+(ya-up)^2),
                           sqrt((xa-left)^2+(ya-down)^2),
                           sqrt((xa-right)^2+(ya-down)^2))/delta)
        TAmax <- trunc(max(sqrt((xa-left)^2+(ya-up)^2),
                           sqrt((xa-right)^2+(ya-up)^2),
                           sqrt((xa-left)^2+(ya-down)^2),
                           sqrt((xa-right)^2+(ya-down)^2))/delta)
        for (TA in TAmin:TAmax) {
          p <- 1
          areaParcial <- areaTotal
          areaTotal <- 0
          for (n in seq(down, up-1, p)) {
            x1 <- 0
            x2 <- 0
            raiz <- 2*n*ya - n^2 - ya^2 + ( delta * (TA + 1) )^2
            if (raiz >= 0) {
              
              s1 <- xa - sqrt(raiz)
              s2 <- xa + sqrt(raiz)
              
              if (s1 < right & s2 > left) {
                
                x1 <- max(s1, left)
                x2 <- min(s2, right)
                areaTotal <- areaTotal + (x2 - x1) * p
              }
            }
          }
          area   <- areaTotal - areaParcial
          temp.dt <- data.table(
            antennaID = stringr::str_pad(as.character(k), 2, pad = '0'), 
            centroidCoord_x = (0.5 + (i - 1)) * tileDim_x, 
            centroidCoord_y = (0.5 + (j - 1)) * tileDim_y, 
            TA = as.integer(TA), area = area)
          output.dt <- rbindlist(list(output.dt, temp.dt))
        } # endFor TA
      } # endFor j
    } # endFor i
    cat('ok.\n')
  }

  return(output.dt)
}
