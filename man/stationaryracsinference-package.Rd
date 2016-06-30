\name{stationaryracsinference-package}
\alias{stationaryracsinference-package}
\alias{stationaryracsinference}
\docType{package}
\title{
\packageTitle{stationaryracsinference}
}
\description{
\packageDescription{stationaryracsinference}
}
\details{
The aim of this package is to include functions for estimating various parameters/properties of RACS from single 2D binary digital images. 
 The functions assume that the image is a digitisation of a stationary and mixing RACS.

 This is particularly relevant to remotely sensed thematic maps and landscape pattern indices (see Hingee 2016 for a discussion).

The package expects RACS observations to be in binary raster \code{owin} format.
 \code{owin} objects are used because \code{spatstat} has many relevant functions for manipulating such objects. 
Conversion of binary raster images in other formats to \code{owin} is quite easy by first converting any images to \code{spatstat}'s \code{im} format (using \code{as.im()}).
 \code{as.owin()} can then be applied to get an owin object.

 Typically the functions require passing both an \code{owin} object containing the random set observation (an image) and an \code{owin} object specifying the boundary of the observation.

 An example application to data in shapefile and ERMapper format is contained in an attached vignette.

The main functions in this (very small) package are:

\code{coveragefrac()} for estimating coverage probability 

\code{contagSphCont()} for estimating disc-state contagion

Another function \code{Hest()} from \code{spatstat} is used heavily in Hingee (2016).

\textbf{Notation for manual:}
Though out this manual the RACS process is referred to as `Xi' or $\Xi$, with the observation (the image) called `xi' or $\xi$.


}
\author{
\packageAuthor{stationaryracsinference}

Maintainer: \packageMaintainer{stationaryracsinference}
}
\references{
Hingee, K.L. (2016) Statistics for Patch Observations. ISPRS Congress Proceedings, IPSRS.
}
\keyword{ package }
\keyword{ spatial }
\examples{
%coverage probability, core probability, SCD, and disc-state contagion estimation
XiOWIN <- heather$coarse
windowOWIN <- Frame(heather$coarse)
 
coverageprob <- coveragefrac(XiOWIN,windowOWIN)

scdest <- Hest(XiOWIN,W=windowOWIN)
coreprob <- Hest(complement.owin(XiOWIN),W=windowOWIN)
discstatecontag <- contagSphCont(scdest,coreprob,coverageprob)
par(mfrow=c(1,3))
plot(scdest)
plot(coreprob)
plot(discstatecontag)
}
