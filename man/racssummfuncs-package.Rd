\name{An overview of this package}
\alias{racssummfuncs-package}
\alias{racssummfuncs}
\docType{package}
\title{
\packageTitle{racssummfuncs}: An Overview
}
\description{
\packageDescription{racssummfuncs}
}
\details{
This package contains functions that estimate various parameters/properties of random closed sets (RACS) from single 2D bi-tonal digital images.
 It enables description of spatial patterns from a rigorous RACS foundation.
 The theoretical background and motivation for these tools can be found in \cite{Hingee (2016)}, \cite{Hingee (2017)} and \cite{Hingee (2018PhDThesis**)}.
 An introduction to using these tools with remotely sensed data is included in  the vignette \emph{remotesensedatademo}. 
 This can be accessed by starting \R and typing \kbd{library(racssummfuncs)} and \kbd{vignette('remotesensedatademo')}.
 Specific help for each function can be accessed by \kbd{help(fun)} where \kbd{fun} is replaced by the name of the function.

The functions accept observations in raster \code{owin} format provided by the \pkg{spatstat} package, or for \command{mvlc} and \code{mvlgb} the observation can be provided as an \code{im} object (also provided by \pkg{spatstat}) with pixel values of only \eqn{1}s, \eqn{0}s and \code{NA}s.
 These formats were chosen because \pkg{spatstat} has many tools that ease the statistical analysis of spatial data, and many of the functions in \pkg{racssummfuncs} take advantage of these. 
 Conversion from other raster formats to \code{im} is quite easy using either \code{\link[spatstat]{im}} or \code{\link[spatstat]{as.im}}.
 The function \code{\link[spatstat]{as.owin}} can then be used to obtain an \code{im} object; note that for the \code{owin} object to correspond to all the \eqn{1}s in the \code{im} object, every non-\eqn{1} pixel must be first set to \code{NA}. 
 
Typically the functions will also require passing an \code{owin} object that specifies the observation window.

The summary functions that this package can generate include:
**include list, function and history/use case.

Additionally \code{\link[spatstat]{Hest}} from \pkg{spatstat} will estimate spherical contact distributions (used heavily in \cite{Hingee (2016)}) which can also be used to estimate the \emph{Minkowski-Bouligand Fractal Dimension} \cite{**still looking**}.

**what about background warning/theory - like the RACS being stationary and mixing**
}
\section{Notation}{
Through out this manual an observation (a bi-tonal \code{im} object, or \code{owin} object with an observation region) is referred to as \sQuote{xi}.
 The underlying RACS from which this observation comes is referred to as \sQuote{Xi}.
}

\author{
\packageAuthor{racssummfuncs}

Maintainer: \packageMaintainer{racssummfuncs}
}
\references{
Hingee, K.L. (2016) Statistics for Patch Observations. ISPRS - International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences pp. 235-242. ISPRS.

Hingee, K.L. (2017) A covariance-based estimator and relationship for lacunarity. In preparation.

Hingee, K.L. (2018) Spatial statistics of random closed sets for Earth observations. The University of Western Australia. PhD Thesis. In preparation.
}
\keyword{ package }
\keyword{ spatial }
\concept{ landscape ecology }
\concept{ pattern metric }
\concept{ landscape pattern indices }
\examples{
%coverage probability, core probability, SCD, and disc-state contagion estimation
xi <- heather$coarse    #an observation of Xi
obswindow <- Frame(heather$coarse)
 
#estimate coverage probability
coverageprob <- coveragefrac(xi,obswindow)

#estimate covariance from xi
xi.covariance <- covariance(xi)
#estimate spherical contact distribution, core probability and disc-state contagion
scdest <- Hest(xi,W=obswindow)   #Hest is from spatstat package
coreprob <- Hest(complement.owin(xi),W=obswindow)
discstatecontag <- contagSphCont(scdest,coreprob,coverageprob)

plot(scdest, main="Conditional Spherical Contact Distribution Estimates\n
 Showing Estimates Using Various Border Corrections")
plot(coreprob, main = "Conditional Core Probability Estimates\n
 Showing Estimates using various Border Corrections")
plot(discstatecontag, main = "A Disc-State Contagion Estimate\n
 (Using the Reduced-Sample Border Correction)")

#Mass-Variance Lacunarity Estimates. If RcppRoll is installed than the gliding box estimate can also be generated for comparison.
xi.mvlc <- mvlc(c(1,2,3,4), covariance = xi.covariance, p = coverageprob)
xi.im <- as.im(xi, na.replace = 0)
xi.mvlgb <- mvlgb(xi.im, c(1,2,3,4))

#Simulation of Boolean Models
booldisc.rfixed.xi <- rBooleanDetermDiscs(2.2064E-3,2,obswindow)
}
