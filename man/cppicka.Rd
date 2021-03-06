% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cppicka.R
\name{cppicka}
\alias{cppicka}
\title{Picka's Reduced Window Estimator of Coverage Probability}
\usage{
cppicka(xi, obswin = NULL, setcov_boundarythresh = NULL)
}
\arguments{
\item{xi}{An observation of a RACS of interest as a full binary map (as an \code{im} object) or as the foreground set (as an \code{owin} object).
In the latter case the observation window, \code{obswin}, must be supplied.}

\item{obswin}{If \code{xi} is an \code{owin} object then \code{obswin} is an
\code{owin} object that specifies the observation window.}

\item{setcov_boundarythresh}{To avoid instabilities caused by dividing by very small quantities, if the set covariance of the observation window
is smaller than \code{setcov_boundarythresh}, then the returned pixel value is NA.}
}
\value{
An \code{im} object. Pixel values correspond to estimates of the coverage probability
from the subregion of the observation window, \eqn{W}, that is the intersection of \eqn{W} and \eqn{W} shifted by vector \eqn{v}, where \eqn{v} is the pixel location.
}
\description{
This function provides estimates of coverage probability from subsets of the observation window,
which are a key component of balanced estimators of covariance, centred covariance, pair-correlation and gliding box lacunarity.
}
\details{
The plug-in moment covariance estimator (\code{\link{plugincvc}}) uses less of the observation window than the usual coverage probability estimators.
Picka (1997, 2000) created new 'balanced' estimators of centred covariance and pair-correlation
that accounted for this difference.
A key component of Picka's estimators is an estimate of the coverage probability from the subregion of the binary map that is
the intersection between \eqn{W} and \eqn{W} shifted by vector \eqn{v}, where \eqn{W} is the observation window (p.~687, Picka, 2000).
If we treat \eqn{X} and \eqn{W} as indicator functions representing the foreground and observation window respectively,
this coverage probability estimator used by Picka is
\deqn{ \frac{\int X(u) W(u) W(u - v) du} {\int W(u) W(u - v) du}. }{integral(X(u) W(u) W(u - v) du)  /  integral(W(u) W(u - v) du).}

\code{cppicka} produces these estimates for an array of vectors \eqn{v} using fast Fourier transforms.
}
\examples{
xi <- heather$coarse
obswindow <- Frame(heather$coarse)
cp <- coverageprob(xi, obswindow)
cpp1 <- cppicka(xi, obswindow)
}
\references{
Picka, J.D. (1997) \emph{Variance-Reducing Modifications for Estimators of Dependence in Random Sets}. Ph.D.: Illinois, USA: The University of Chicago.

Picka, J.D. (2000) Variance reducing modifications for estimators of standardized moments of random sets. \emph{Advances in Applied Probability}, 32, 682-700.
}
\author{
Kassel Liam Hingee
}
\keyword{nonparametric}
\keyword{spatial}
