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
The aim of this package is to include functions for estimating various parameters/properties of random closed set (RACS) from single 2D binary digital images. 
 This is particularly relevant to remotely sensed thematic maps and landscape pattern indices.
 More background and details can be found in Hingee (2016).

 The estimators represent the image as a digitisation of a stationary and mixing RACS.
 A RACS is any process that generates closed sets (intuitively one may think of them as generating sets that include their boundary).
 A stationary and mixing RACS is a RACS for which
 the statistical properties are invariant to translation (stationarity) and observations at locations far away from each other approach independence (mixing). 

Thus this package provides a way to describe a statistical representation of a landscape. 
These are more abstract and synoptic descriptions (and thus hopefully less sensitive to measurement artefacts) than Landscape Pattern Indices.


The package expects observations to be in binary raster \code{owin} format (this format is provided by the \code{spatstat} package).
 \code{owin} objects are used because \code{spatstat} has many relevant tools for statistical analysis of spatial data. 
Conversion of binary raster images in other formats to \code{owin} is quite easy by first converting any images to \code{spatstat}'s \code{im} format (using \code{as.im()}).
 \code{as.owin()} can then be applied to get an owin object.
 Typically the functions require passing both an \code{owin} object containing the random set observation (an image) and an \code{owin} object specifying the boundary of the observation.

 An example application to data in shapefile and ERMapper format is contained in an attached vignette \emph{remotesensedatademo}. 
 This can be accessed by starting \R and typing \code{library(stationaryracsinference)} and \code{vignette('remotesensedatademo')}.

The only functions in this (very small) package are:
\itemize{
\item \code{coveragefrac()} for estimating coverage probability 
\item \code{contagSphCont()} for estimating disc-state contagion
}

Another function \code{Hest()} from \code{spatstat} is used heavily in Hingee (2016) to estimate the spherical contact distributions and conditional core probability.

}
\section{Notation}{
Through out this manual the RACS process is referred to as `Xi' or \eqn{\Xi}, with the observation (the image) called `xi' or \eqn{\xi}.
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
xi <- heather$coarse
obswindow <- Frame(heather$coarse)
 
coverageprob <- coveragefrac(xi,obswindow)

scdest <- Hest(xi,W=obswindow)
coreprob <- Hest(complement.owin(xi),W=obswindow)
discstatecontag <- contagSphCont(scdest,coreprob,coverageprob)

plot(scdest, main="Conditional Spherical Contact Distribution Estimates\n Showing Estimates Using Various Border Corrections")
plot(coreprob, main = "Conditional Core Probability Estimates\n Showing Estimates using various Border Corrections")
plot(discstatecontag, main = "A Disc-State Contagion Estimate\n (Using the Reduced-Sample Border Correction)")
}
