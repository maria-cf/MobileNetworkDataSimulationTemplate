#' @title Compute equivalence between tiles and raster cells.
#'
#' @description This function computes the equivalence in the index of each 
#' raster cell (according to the raster package) and each tile (according to the
#' network event data simulator). Raster cells are counted from left to right 
#' starting in the upper left corner by 1. Tiles are counted from left to right 
#' starting in the bottom left corner by 0.
#'
#' @param nrow Number of rows in the raster.
#' 
#' @param ncol Number of columns in the raster.
#' 
#' @param order_by Either \code{rasterCellID} (default) or \code{tileID}.
#' 
#' @return A matrix with two named columns (\code{rasterCellID} and 
#' \code{tileID}) with the equivalence for each raster cell/tile.
#'
#' @examples
#' tileEquivalence(4, 3)
#' 
#' tileEquivalence(4, 3, order_by = 'tileID')  
#'
#' @export
tileEquivalence <- function(nrow, ncol, order_by = c('rasterCell')){

  nCells <- nrow * ncol
  tileID_raster <- as.vector(matrix(1:nCells, ncol = ncol, byrow = TRUE))
  tileID_simulator <- as.vector(matrix(0:(nCells - 1), ncol = ncol, byrow = TRUE)[nrow:1, ])
  tileCorresp <- cbind(rasterCell = tileID_raster, tile = tileID_simulator)
  tileCorresp <- tileCorresp[order(tileCorresp[, order_by]), ]
  return(tileCorresp)

}
