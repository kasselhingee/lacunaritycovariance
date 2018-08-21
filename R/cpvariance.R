#' @title Variances and Confidence Intervals for Coverage Fraction Estimates 
#' @aliases asympvarP
#' @export cpvariance varCovProb_ests cpvariance.covarsupplied
#' 
#' @description Functions for estimating variance (and confidence intervals)  of the coverage fraction estimates
 
#' @details Molchanov (1997 Ch3.1) notes that the variance of the estimator \code{\link[stationaryracsinference]{coveragefrac}} is 
#' \deqn{var = \frac{1}{A(w)^2}\int_W \gamma_W(v)C(v) - p^2 dv}
 #estimate exact variance of pest
 #from set covariance of W is the erosion of W by {o,v}
 #exact var = 1/A(Window)^2*sumOverVectorsvOF[setcov(Xi,v)-setcov(W,v)*p^2]


#' @details This function assumes that the probability of two points being in \eqn{\Xi} is independent \deqn{C(v) = p^2} outside the calculated window of covariance which is...


#' @return varcovProb returns a variance estimate using the covariance function and an estimate of area fraction). This will not be as good as the spectral density based estimate**.
#' asympvarP() isn't constructed yet but will estimate a variance assuming that things are close to Gaussian distributions.
#' @param Xi is an observation (in owin) format of a RACS or an image of 0, 1 and NA values.
#' @param obswin is the corresponding observation window. 
#' @param covar A covariance image for the RACS (could be estimated from Xi)
#' @param setcov_boundarythresh When estimating covariance the threshold at which the set covariance of the observation window is deemed too small and estimation doesn't occur. See \code{tradcovarest()}.
#' @param modifications A list of modifications of centred covariance estimation to use - see \code{ccvc()}.

#' @examples 
#' Xi <- heather$coarse
#' obswin <- Frame(Xi)
#' varCovProb_ests(Xi, obswin, modifications = "all")
#' @references 
#' Molchanov, I. (1997) Statistics of the Boolean Model for Practitioners and Mathematicians. Wiley.
cpvariance <- function(Xi, obswin){
   Xiinside <- intersect.owin(Xi,obswin)
   setcovXi <- setcov(Xiinside)
   setcovB <- setcov(obswin)
   p <- max(setcovXi)/max(setcov(obswin)) #using this instead of normal phat estimate seems to make positive values more likely at least? **I'd really like to know why!
   integrand <- eval.im(setcovXi-(p^2)*setcovB, harmonize = TRUE)
   #test that integrand reaches 0
   edgeValues = c(integrand[1,-1],integrand[-1,ncol(integrand)],
                  integrand[nrow(integrand),-ncol(integrand)],integrand[-nrow(integrand),1])
   if (max(edgeValues,na.rm=TRUE) > diff(range(integrand,na.rm=TRUE))*1e-5){
     warning("covariance weighted by set covariance of the window isn't uniformly close to p^2 at boundary\n")
     cat("max size of (C(v) - p^2)*[Set Covariance of Window] on boundary is ", max(edgeValues,na.rm=TRUE),"\n",sep="")
   }
   return((1/(area.owin(obswin))^2)*sum(integrand)*integrand$xstep*integrand$ystep)
 }



#a seperate function could be useful because the othe function will have less machine error
#' @describeIn cpvariance Variance estimate from a given covariance function
cpvariance.covarsupplied <- function(covar, obswin){
  p <- covar[as.ppp(c(0,0), W = Frame(covar))]
  setcovB <- setcov(obswin)
  integrand <- eval.im((covar-(p^2))*setcovB, harmonize = TRUE)
  #test that integrand reaches 0
  edgeValues = c(integrand[1,-1],integrand[-1,ncol(integrand)],
                integrand[nrow(integrand),-ncol(integrand)],integrand[-nrow(integrand),1])
  if (max(edgeValues,na.rm=TRUE) > diff(range(integrand,na.rm=TRUE))*1e-5){
   warning("covariance weighted by set covariance of the window isn't uniformly close to p^2 at boundary\n")
   cat("max size of (C(v) - p^2)*[Set Covariance of Window] on boundary is ", max(edgeValues,na.rm=TRUE),"\n",sep="")
  }
  return((1/(area.owin(obswin))^2)*sum(integrand)*integrand$xstep*integrand$ystep)
} 


 #from molchanov (and the limit of the above) as A(W)--> infy, var(p)-->0 which makes sense
 #however if we look psqrt(A(W)) then something different happens?? *ASK GOPAL
 asympvarP <- function(){}

#' @describeIn cpvariance Use multiple balanced estimators of covariance to estimate variance of coverage probability
varCovProb_ests <- function(Xi, obswin = NULL,
        setcov_boundarythresh = NULL,
        modifications = "all"){
  cvchat <- tradcovarest(Xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  cpp1 <- cppicka(Xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  phat <- cvchat[ppp(x = 0, y = 0, window = Frame(cvchat))]
  #phat <- coverageprob(Xi, obswin)
  
  cvchats <- racscovariance.cvchat(cvchat, cpp1, phat, modifications = modifications, drop = FALSE) 
  
  if (is.null(obswin) && is.im(Xi)){
    obswin <- as.owin(Xi) #only excludes NA values in Xi
  }
  if (is.mask(Xi) || is.im(Xi)){  setcovW <- setcov(obswin, xy = Xi) 
  } else { setcovW <- setcov(obswin, xy = cvchat) }
  setcovW <- as.im(setcovW, xy = cvchat) #harmonise results

  integrands <- solapply(cvchats, function(x) eval.im(setcovW * (x - phat^2)))
  varests <- vapply(integrands, integral.im, 0.13)/(area.owin(obswin)^2)
   
  return(varests)
}


cpvariance.trad <- function(Xi, obswin){
   Xiinside <- intersect.owin(Xi,obswin)
   setcovXi <- setcov(Xiinside)
   setcovB <- setcov(obswin)
   p <- max(setcovXi)/max(setcov(obswin)) #using this instead of normal phat estimate seems to make positive values more likely at least? **I'd really like to know why!
   integrand <- eval.im(setcovXi-(p^2)*setcovB, harmonize = TRUE)
   #test that integrand reaches 0
   edgeValues = c(integrand[1,-1],integrand[-1,ncol(integrand)],
                  integrand[nrow(integrand),-ncol(integrand)],integrand[-nrow(integrand),1])
   if (max(edgeValues,na.rm=TRUE) > diff(range(integrand,na.rm=TRUE))*1e-5){
     warning("covariance weighted by set covariance of the window isn't uniformly close to p^2 at boundary\n")
     cat("max size of (C(v) - p^2)*[Set Covariance of Window] on boundary is ", max(edgeValues,na.rm=TRUE),"\n",sep="")
   }
   return((1/(area.owin(obswin))^2)*sum(integrand)*integrand$xstep*integrand$ystep)
 }
