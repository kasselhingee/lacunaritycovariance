
#Estimate the coverage probability of Xi by calculating the area of Xi (in pixels) within a given boundary and dividing by the area
#in reality the area of Xi is probably slightly different - depends on the sensing method of Xi. Can assume really close though.
 covpest <- function(Xi,w){
   stopifnot(is.owin(Xi))
   stopifnot(is.owin(w))   
   XiInsideW <- intersect.owin(Xi,w)
   areaXiInside <- area.owin(XiInsideW)
   areaWindow <- area.owin(w)
   
   covProbEstimate <- areaXiInside/areaWindow
   return(covProbEstimate)
 }
 
#estimate exact variance of pest
#from Molchanov Ch3.1 and that the set covariance of W is the erosion of W by {o,v}
#exact var = 1/A(Window)*sumOverVectorsVOF[setcovariance(W,v)*(covarianceEst(v) - p^2)]
#assume that covariance = p^2 outside the covariance map. BUT give an error if covariance Map is too large at edges?
#assume covarianceMap has different slots
 varPexactEst <- function(pest,covarianceMap,w){
   ##need some checking that covariance=p^2 outside covariance map
   edgeValues = c(covarianceMap$covariance[1,-1],covarianceMap$covariance[-1,ncol(covarianceMap$covariance)],
                  covarianceMap$covariance[nrow(covarianceMap$covariance),-1],covarianceMap$covariance[-1,1])
   if (max(abs(edgeValues-pest^2),na.rm=TRUE) > diff(range(covarianceMap$covariance,na.rm=TRUE))*1e-5){
     warning("covariance map doesn't reached p^2 at boundary")
   }
   return(1/(area.owin(w)^2)*vectorSteps*sum(covarianceMap$numerator-covarianceMap$denominator*pest^2,na.rm=TRUE))
 }
 

