% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbdr.R
\name{rbdr}
\alias{rbdr}
\alias{bdrcoverageprob}
\alias{bdrcovar}
\title{Simulation of Boolean Model of Deterministic Rectangles}
\usage{
rbdr(lambda, grain, win, seed = NULL)

bdrcoverageprob(lambda, grain)

bdrcovar(lambda, grain, xy)
}
\arguments{
\item{lambda}{Intensity of the germ process (which is a Poisson point process)}

\item{grain}{Rectangle object specifying the grain}

\item{win}{The window to simulate in (an \code{owin} object)}

\item{seed}{Optional input (default in NULL). Is an integer passed to \code{\link[base]{set.seed}}. Used to reproduce patterns exactly.}

\item{xy}{A raster object that specifies the pixel coordinates of the desired covariance image. \code{xy} works in similar fashion to passing an image or pixel mask through the \code{xy} argument of \code{\link[spatstat.geom]{as.mask}} in \pkg{spatstat}.}
}
\value{
Depends on the function used (see Functions section).
}
\description{
Functions for simulating a Boolean model with grains that are deterministic rectangles.
A Boolean model is a two stage model, first the locations (called germs) of grains are randomly distributed according to a Poisson point process, then a random grain is placed on each germ independently.
An introduction can be found in (Chiu et al., 2013).
Also described in this help file are functions for calculating the coverage probability and covariance.
}
\section{Functions}{
\itemize{
\item \code{rbdr()}: Returns an \code{owin} that is a set generated by simulating a Boolean
model with a specified intensity and fixed rectangular grain.
The window information is not contained in this object.
If the simulated set is empty then an empty \code{owin} object is returned.
The point process of germs is generated using spatstat's \code{\link[spatstat.random]{rpoispp}}.

\item \code{bdrcoverageprob()}: Returns the true coverage probability given the intensity and grain.

\item \code{bdrcovar()}: Returns an image of the covariance as calculated from disc radius and intensity.

}}
\section{WARNING}{

The returned object of \code{rbdr} is only the foreground of a binary map and thus can have much smaller extent than the simulation window (e.g. when the simulated set is empty).
}

\examples{
grain <- owin(xrange = c(-5, 5), yrange = c(-5, 5))
win <- owin(xrange = c(0, 100), c(0, 100))
lambda <- 4.2064E-3
xi <- rbdr(lambda, grain, win)

cp_theoretical <- bdrcoverageprob(lambda, grain)
xy <- as.mask(dilationAny(win, win), eps = c(1, 1))
cvc_theoretical <- bdrcovar(lambda, grain, xy)
}
\references{
Chiu, S.N., Stoyan, D., Kendall, W.S. and Mecke, J. (2013) \emph{Stochastic Geometry and Its Applications}, 3rd ed. Chichester, United Kingdom: John Wiley & Sons.
}
\keyword{datagen}
\keyword{spatial}
