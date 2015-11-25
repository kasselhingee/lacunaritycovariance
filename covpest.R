
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
 #exact var = 1/A(Window)^2*sumOverVectorsvOF[setcovariance(W,v)*(covarianceEst(Xi,v) - p^2)]
 #          = 1/A(Window)^2*sumOverVectorsvOF[setcov(Xi,v)-setcov(W,v)*p^2]
 #assume that covariance = p^2 outside the covariance map. BUT give an error if covariance Map is too large at edges?
 exactvarP <- function(XiOWIN,w){
   Xiinside <- intersect.owin(XiOWIN,w)
   p <- area.owin(Xiinside)/area.owin(w)
   setcovXi <- setcov(Xiinside)
   setcovB <- setcov(w)
   integrand <- eval.im(setcovXi-p^2*setcovB,harmonize = TRUE)
   #test that integrand reaches 0
   edgeValues = c(integrand[1,-1],integrand[-1,ncol(integrand)],
                  integrand[nrow(integrand),-ncol(integrand)],integrand[-nrow(integrand),1])
   if (max(integrand,na.rm=TRUE) > diff(range(integrand,na.rm=TRUE))*1e-5){
     warning("covariance isn't uniformly close to p^2 at boundary\n")
     cat("max difference between C(v) and p^2 on boundary is ", max(integrand,na.rm=TRUE),sep="")
   }
   return((1/(area.owin(w))^2)*sum(integrand)*integrand$xstep*integrand$ystep)
 }
 
 
 
#estimate exact variance of pest
#from Molchanov Ch3.1 and that the set covariance of W is the erosion of W by {o,v}
#exact var = 1/A(Window)*sumOverVectorsVOF[setcovariance(W,v)*(covarianceEst(v) - p^2)]
#assume that covariance = p^2 outside the covariance map. BUT give an error if covariance Map is too large at edges?
 #assume covarianceMap has different slots
 exactvarP_direct <- function(pest,covarianceMap,w){
   ##need some checking that covariance=p^2 outside covariance map
   edgeValues = c(covarianceMap$covariance[1,-1],covarianceMap$covariance[-1,ncol(covarianceMap$covariance)],
                  covarianceMap$covariance[nrow(covarianceMap$covariance),-ncol(covarianceMap$covariance)],covarianceMap$covariance[-nrow(covarianceMap$covariance),1])
   if (max(abs(edgeValues-pest^2),na.rm=TRUE) > diff(range(covarianceMap$covariance,na.rm=TRUE))*1e-5){
     warning("covariance map doesn't reached p^2 at boundary")
   }
   return(1/(area.owin(w)^2)*vectorSteps*sum(covarianceMap$numerator-covarianceMap$denominator*pest^2,na.rm=TRUE))
 }
 
 
 
 
 #from molchanov (and the limit of the above) as A(W)--> infy, var(p)-->0 which makes sense
 #however if we look psqrt(A(W)) then something different happens?? *ASK GOPAL
 asympvarP <- function(){}

 #calculates width of symmetric confidence interval given a required confidence level and a sd
 confidenceIntervalOfGaussian <- function(confidenceLevel,sd){
   minBound <- qnorm((1-confidenceLevel)/2,sd=sd)
   maxBound <- -minBound
   return(c(minBound,maxBound))
 }