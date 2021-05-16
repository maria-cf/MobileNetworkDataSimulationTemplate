#' @title Read jointly all output files from the network data event simulator.
#'
#' @description This function reads all output files from the network data event 
#' simulator from an output path name specified by the user and returns an rds 
#' file in another input path name for each of the files.
#'
#' @param ntiles_x number of tiles in axis X.
#'
#' @param ntiles_y number of tiles in axis Y.
#' 
#' @param tile_sizeX the size in meters of the tile in axis X.
#' 
#' @param tile_sizeY the size in meters of the tile in axis Y.
#' 
#' @param tileEquiv.dt data.table with the equivalence between tiles and raster cells.
#' 
#' @return A data.table with columns: 
#'
#' @examples
#' \dontrun{
#' 
#' }
#'
#' @import data.table
#'
#' @include tileEquivalence.R
#'
#' @export
build_centroid_regions <- function(ntiles_x, ntiles_y, tile_sizeX, tile_sizeY,
                                   tileEquiv.dt){
  
  # # TechDebt: generalize with breaks_x and breaks_y in the input arguments.
  # breaks_x <- c(3000, 5000, 7500)
  # breaks_y <- c(3500)
  
  ntiles     <- ntiles_x * ntiles_y
 
   #### * Build centroid distance matrix                                     ####
  cat('     centroid-distance matrix...')
  centroidCoord_x  <- (0.5 + 0:(ntiles_x - 1)) * tile_sizeX
  centroidCoord_y  <- (0.5 + 0:(ntiles_y - 1)) * tile_sizeY
  centroidCoord.dt <- as.data.table(expand_grid(centroidCoord_y, centroidCoord_x))[
    , 2:1][
      , tile := 0:(ntiles-1)][tileEquiv.dt, on = 'tile']
  cat(' ok.\n')
  
  #### * Build regions, zones, macrozones, total                            ####
  cat('     regions, zones, macrozones, total...')
  centroidCoord.dt[centroidCoord_x <= 3000 & centroidCoord_y <= 3500, region := 1]
  centroidCoord.dt[centroidCoord_x <= 3000 & centroidCoord_y >  3500 & centroidCoord_y <= 6000, region := 2]
  centroidCoord.dt[centroidCoord_x <= 3000                           & centroidCoord_y >  6000, region := 3]
  centroidCoord.dt[centroidCoord_x  > 3000 & centroidCoord_x <= 5000 & centroidCoord_y < 5000, region := 4]
  centroidCoord.dt[centroidCoord_x  > 3000 & centroidCoord_x <= 5000 & centroidCoord_y >= 5000, region := 5]
  centroidCoord.dt[centroidCoord_x  > 5000 & centroidCoord_x <= 7500 & centroidCoord_y < 5000, region := 6]
  centroidCoord.dt[centroidCoord_x  > 5000 & centroidCoord_x <= 7500 & centroidCoord_y >= 5000, region := 7]
  centroidCoord.dt[centroidCoord_x  > 7500 & centroidCoord_y <  2000, region := 8]
  centroidCoord.dt[centroidCoord_x  > 7500 & centroidCoord_y >= 2000 & centroidCoord_y < 7000, region := 9]
  centroidCoord.dt[centroidCoord_x  > 7500 & centroidCoord_y >= 2000 & centroidCoord_y >= 7000, region := 10]
  
  centroidCoord.dt[region %in% c(1, 2, 3),  zone := 1]
  centroidCoord.dt[region %in% c(4, 6),     zone := 2]
  centroidCoord.dt[region %in% c(5, 7),     zone := 3]
  centroidCoord.dt[region %in% c(8, 9, 10), zone := 4]
  
  centroidCoord.dt[zone %in% c(1, 2),  macrozone := 1]
  centroidCoord.dt[zone %in% c(3, 4),  macrozone := 2]
  
  centroidCoord.dt[macrozone %in% c(1, 2),  total := 1]
  cat(' ok.\n')
  
  return(centroidCoord.dt)
  
  
}
  