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


The only functions in this (very small) package are:
\itemize{
\item \code{coveragefrac()} for estimating coverage probability 
\item \code{contagSphCont()} for estimating disc-state contagion
}

Another function \code{Hest()} from \code{spatstat} is used heavily in Hingee (2016) to estimate the spherical contact distributions and conditional core probability.

}
\section{Notation}{
Through out this manual the RACS is referred to as `Xi' and an observation of Xi (which is a bitonal image) is referred to as `xi'.
}

\author{
\packageAuthor{racssummfuncs}

Maintainer: \packageMaintainer{racssummfuncs}
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

plot(scdest, main="Conditional Spherical Contact Distribution Estimates\n
 Showing Estimates Using Various Border Corrections")
plot(coreprob, main = "Conditional Core Probability Estimates\n
 Showing Estimates using various Border Corrections")
plot(discstatecontag, main = "A Disc-State Contagion Estimate\n
 (Using the Reduced-Sample Border Correction)")
}
