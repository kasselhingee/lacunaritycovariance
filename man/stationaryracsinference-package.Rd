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
 In more rigorous language the functions in this package assume that the image is a digitisation of a stationary and mixing RACS. 
 These functions estimate various properties of RACS from the image.

 This is particularly relevant to remotely sensed thematic maps and landscape pattern indices (Hingee 2016).

The package expects RACS observations to be in binary raster \code{owin} format.
 \code{owin} objects are used because spatstat has many relevant functions for manipulating such objects. 
 Typically functions require passing both an owin object containing the random set and an owin object specifying the boundary of the observation.

Conversion of binary raster images in other formats to \code{owin} is quite easy
by first converting any images to \code{spatstat}'s \code{im} format (using \code{as.im()}).
 \code{as.owin()} can then be applied to get an owin object.


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
Hingee, K.L. (2016) Statistics for Patch Observations. ISPRS Congress Proceedings p. IPSRS.

refer to Chiu et al 2013, and Molchanov**
}
\keyword{ package }
\keyword{ spatial }
\examples{
%coverage and covariance estimation
XiOWIN <- heather$coarse
windowOWIN <- Frame(heather$coarse)
 
coverageProb <- coveragefrac(XiOWIN,windowOWIN)
 
covarianceFcn <- covariance(XiOWIN,windowOWIN)
plot(covarianceFcn$covariance)
plot(covarianceFcn$covariance - coverageProb*coverageProb)
}
