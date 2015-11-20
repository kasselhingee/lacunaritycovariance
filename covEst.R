
#estimate the covariance, or two point probability
#inputs must be spatstat owin or list of two elements
#Input: a map of Xi, possible extending outside the boundary
#       a boundary
#       a vector in (x,y) coordinates (any real values work, but I suspect only multiples of the resolution would make sense)
covarianceEst <- function(Xi,boundary,v){
  denominator <- area.owin(intersect.owin(boundary,shift.owin(boundary,vec=v)))#need to handle denominator of 0
  if (denominator == 0){
    covarianceEst <- NA
    return(list(numerator = NA, denominator = NA, covarianceEst = NA))
    }
  else {  
    Xiinside <- intersect.owin(Xi,boundary) 
    numerator <- area.owin(intersect.owin(Xiinside,shift.owin(Xiinside,vec=v)))
    covarianceEst <- numerator/denominator
    }
  
  return(list(numerator = numerator, denominator = denominator, covarianceEst = covarianceEst))
}



#estimate a map of covariance for each vector - pretty naive computationally Xi in OWIN, boundary in OWIN, maxshifts in units the same as Xi 
#ignores point estimates that use an area smaller than 10% of the window
covarianceMapEst <- function(Xi,boundary,maxXshiftdistance,maxYshiftdistance){
  windowArea <- area.owin(boundary)
  #create the vectors for testing
  shiftVectorX <- c(-rev(seq(0,maxXshiftdistance,by=Xi$xstep)),seq(Xi$xstep,maxXshiftdistance,by=Xi$xstep))
  shiftVectorY <- c(-rev(seq(0,maxYshiftdistance,by=Xi$ystep)),seq(Xi$ystep,maxYshiftdistance,by=Xi$ystep))
  covarMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
  numeratorMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
  denominatorMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
  for (i in 1:length(shiftVectorX)){
    for (j in 1:length(shiftVectorY)){
      covarianceEstimation <- covarianceEst(XiOWIN,boundaryOWIN,c(shiftVectorX[i],shiftVectorY[j]))
      if (is.na(covarianceEstimation$denominator) || (covarianceEstimation$denominator < 0.1*windowArea)){
        covarMap[i,j] <- NA
      }
      else {covarMap[i,j] <- covarianceEstimation$covarianceEst}
      numeratorMap[i,j] <- covarianceEstimation$numerator
      denominatorMap[i,j] <- covarianceEstimation$denominator
    }
  } 
  
  covarianceMap = list(covariance = covarMap, numerator=numeratorMap, denominator = denominatorMap, xcol = shiftVectorX, ycol = shiftVectorY, xstep=Xi$xstep,ystep=Xi$ystep)
  return(covarianceMap) 
}