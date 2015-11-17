
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