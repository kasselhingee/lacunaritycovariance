#' @title Variances and Confidence Intervals for Coverage Fraction Estimates 
#' @aliases asympvarP
#' @export varCovProb varCovProb.covarsupplied
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
#' @param Xi is an observation (in owin) format of a RACS.
#' @param w is the corresponding observation window.

#' @examples 
#' #**To come later** 
#' @references 
#' Molchanov, I. (1997) Statistics of the Boolean Model for Practitioners and Mathematicians. Wiley.
 varCovProb <- function(Xi, w){
   Xiinside <- intersect.owin(Xi,w)
   setcovXi <- setcov(Xiinside)
   setcovB <- setcov(w)
   p <- max(setcovXi)/max(setcov(w)) #using this instead of normal phat estimate seems to make positive values more likely at least? **I'd really like to know why!
   integrand <- eval.im(setcovXi-(p^2)*setcovB, harmonize = TRUE)
   #test that integrand reaches 0
   edgeValues = c(integrand[1,-1],integrand[-1,ncol(integrand)],
                  integrand[nrow(integrand),-ncol(integrand)],integrand[-nrow(integrand),1])
   if (max(edgeValues,na.rm=TRUE) > diff(range(integrand,na.rm=TRUE))*1e-5){
     warning("covariance weighted by set covariance of the window isn't uniformly close to p^2 at boundary\n")
     cat("max size of (C(v) - p^2)*[Set Covariance of Window] on boundary is ", max(edgeValues,na.rm=TRUE),"\n",sep="")
   }
   return((1/(area.owin(w))^2)*sum(integrand)*integrand$xstep*integrand$ystep)
 }
  
#a seperate function could be useful because the othe function will have less machine error
#' @describeIn varCovProb Variance estimate from a given covariance function
varCovProb.covarsupplied <- function(covar, w){
  p <- covar[as.ppp(c(0,0), W = Frame(covar))]
  setcovB <- setcov(w)
  integrand <- eval.im((covar-(p^2))*setcovB, harmonize = TRUE)
  #test that integrand reaches 0
  edgeValues = c(integrand[1,-1],integrand[-1,ncol(integrand)],
                integrand[nrow(integrand),-ncol(integrand)],integrand[-nrow(integrand),1])
  if (max(edgeValues,na.rm=TRUE) > diff(range(integrand,na.rm=TRUE))*1e-5){
   warning("covariance weighted by set covariance of the window isn't uniformly close to p^2 at boundary\n")
   cat("max size of (C(v) - p^2)*[Set Covariance of Window] on boundary is ", max(edgeValues,na.rm=TRUE),"\n",sep="")
  }
  return((1/(area.owin(w))^2)*sum(integrand)*integrand$xstep*integrand$ystep)
} 


 #from molchanov (and the limit of the above) as A(W)--> infy, var(p)-->0 which makes sense
 #however if we look psqrt(A(W)) then something different happens?? *ASK GOPAL
 asympvarP <- function(){}



