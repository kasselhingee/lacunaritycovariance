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
This package expects RACS observations to be in the form of binary raster images.
Once a binary image is in R (as an sp.image) then coverage fraction and other things can be calculated.

The aim is to develop methods for estimating various parameters for RACS models of tree canopy and other remotely sensed thematic maps.

The main functions are:

coveragefrac() for coverage fraction

covariance() for the covariance (also known as the two-point probability function)

}
\author{
\packageAuthor{stationaryracsinference}

Maintainer: \packageMaintainer{stationaryracsinference}
}
\references{
refer to Chiu et al 2013, and Molchanov**
}
\keyword{ package }
\keyword{ spatial }
\examples{
#coverage and covariance estimation
XiOWIN <- heather$coarse
windowOWIN <- Frame(heather$coarse)
 
coverageProb <- coveragefrac(XiOWIN,windowOWIN)
 
covar <- balancedracscovariances(XiOWIN,windowOWIN, modification = "pickaH")
}
