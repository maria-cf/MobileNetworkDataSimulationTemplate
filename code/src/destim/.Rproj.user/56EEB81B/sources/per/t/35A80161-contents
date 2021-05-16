#' @title Higher-order level wrapper for HMM fitting functions.
#'
#' @description This is a higher-order wrapper for different functions of 
#' package destim to fit an HMM and compute location probabilities all together.
#' 
#' @param model Object of class \code{HMM} with the initialized HMM model.
#' 
#' @param observedValues Vector with the observed values of the event variables.
#' 
#' @param pad_coef Integer coefficient denoting the number of time periods 
#' needed for the time padding.
#' 
#' @param init logical variable passed to the parameter \code{init} of function
#' \code{fit} to fit the HMM. See ?destim::fit.
#' 
#' @param sstates logical variable indicating whether to compute posterior 
#' location probabilities.
#' 
#' @param scpstates logical variable indicating whether to compute joint 
#' posterior location probabilities.
#'
#' @return A named \code{list} with components model_devID, postLocP, and 
#' postJointLocP with (i) the fitted model (class \code{HMM}), (ii) the 
#' posterior location probabilities (matrix of dimensions tiles $\times$ times),
#'  and (iii) the joint posterior location probabilities (matrix of dimensions 
#'  tiles $\times$ tiles * (times - 1)).
#' 
#' @examples
#' \dontrun{ 
#' }
#'
#' @import data.table
#' 
#' @include fit initparams sstates scpstates
#'
#' @export
compute_HMM <- function(model, observedValues, 
                        pad_coef = 1, init = TRUE, 
                        sstates = TRUE, scpstates = TRUE){
  
  #### Checks over the input arguments                                                          ####
  
  # Check if model is an object of class HMM
  if (!"HMM" %in% class(model)) {
    
    stop(paste0('[compute_HMM] model must be of class HMM.\n'))
  }
  
  emissionProbMatrix <- model$emissions
 
  # Check if emissionProbMatrix meets the restrictions
  if (!all(round(apply(emissionProbMatrix, 1, sum), 6)==1 || round(apply(emissionProbMatrix, 1, sum), 6)==0)) {
    
    stop(paste0('[compute_HMM] There is some problem of consistency in emissionProbMatrix.\n'))
  }
  
  # Check if all observedValues are in emissionProbMatrix
  columnsLack <- setdiff(observedValues, colnames(emissionProbMatrix))
  columnsLack <- columnsLack[!is.na(columnsLack)]
  if (length(columnsLack) > 0) {
    
    stop(paste0('[compute_HMM] The following observedValues are missing in emissionProbMatrix: ', 
                paste0(columnsLack, collapse = ', '),  '.\n'))
    
  }
  
  # Check if is.numeric pad_coef
  if (!is.numeric(pad_coef)) {
    
    stop(paste0('[check_jumps] pad_coef must be numeric.\n'))
  }
  
  ## Time padding
  cat('       time padding...')
  observedValues_pad <- rep(NA, pad_coef * length(observedValues))
  observedValues_pad[seq(1, length(observedValues_pad), by = pad_coef)] <- observedValues
  
  colEvents <- sapply(observedValues_pad, 
                      function(x) ifelse(!is.na(x), 
                                         which (x == colnames(emissionProbMatrix)), NA))
  
  cat(' ok.\n')
  
  #####                             Fit HMM                                #####
  cat('       fitting HMM...')
  fitTry <- try(model_devID <- fit(model, colEvents, init = init)) # ML estimation of transition probabilities
  if(inherits(fitTry, "try-error")){
    # iter <- 1
    # while(inherits(fitTry, "try-error") & iter < 5) {
    #   
    #   model <- initparams(model) # prob trans inicial
    #   model <- initsteady(model) 
    #   fitTry <- try(model_devID <- fit(model, colEvents, init = init)) 
    #   iter <- iter + 1
    #   
    # }
    if(inherits(fitTry, "try-error")){
      stop("Fit model fails")
    }
  }
  cat(' ok.\n')
  
  if (sstates) {
      
    cat('       computing posterior location probabilities...')
    ssTry <- try(A <- sstates(model_devID, colEvents))
      
      if(inherits(ssTry, "try-error")) {
        
        stop("[compute_HMM] Smooth States fails")
      
      }
    cat( ' ok.\n')
    #A <- data.table(as.matrix(A))
    
  } else {
    
    A <- data.table(NULL)
    
  }
  
  if (scpstates) {
    
    cat('       computing joint posterior location probabilities...')
    B <- scpstates(model_devID, colEvents)
    # B <- as.matrix(data.table(B))
    cat(' ok.\n')
    
  } else {
    
    # B <- data.table(NULL)
    
  } 
  output <- list(model_devID = model_devID, postLocP = A, postJointLocP = B)
  
  rm(model_devID)
  rm(A)
  rm(B)
  gc()
  
  return(output)
    
}