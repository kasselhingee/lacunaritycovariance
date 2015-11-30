\name{covarianceEstAtPoint}
\alias{covarianceEstAtPoint}
\alias{covarianceMapEst_direct}
\alias{covarianceRACS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Covariance Estimates - probability of two points being covered%%  ~~function to do ... ~~
}
\description{
These functions estimate the covariance of a stationary random closed set. 
The covariance is also known as the two-point coverage probability, and very closely related to the semivariogram.
 The covariance of a vector $v$ is the probability of two points separated by a $v$ being covered by $\Xi$.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
covarianceEstAtPoint(Xi, boundary, v)
covarianceMapEst_direct(Xi, boundary, maxXshiftdistance, maxYshiftdistance)
covarianceRACS(Xi, boundary, setCovBoundaryThresh = 0.1 * area.owin(boundary))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xi}{
  An observations (in spatstat OWIN format) of the RACS of interest. 
%%     ~~Describe \code{Xi} here~~
}
  \item{boundary}{
  The boundary of the observation in OWIN format. This is needed for the cases where the observation is not a rectangular region.
%%     ~~Describe \code{boundary} here~~
}
  \item{v}{
  A 2D vector in c(x,y) format.
%%     ~~Describe \code{v} here~~
}
  \item{maxXshiftdistance}{The maximum size of x-component of vectors $v$ to estimate
%%     ~~Describe \code{maxXshiftdistance} here~~
}
  \item{maxYshiftdistance}{The maximum size of y-component of vectors $v$ to estimate
%%     ~~Describe \code{maxYshiftdistance} here~~
}
  \item{setCovBoundaryThresh}{To avoid instabilities of dividing by very small areas, any vector $v$ set covariance of the boundary smaller than this threshold is given a covariance of NA}
%%     ~~Describe \code{setCovBoundaryThresh} here~~
}
}
\details{
\code{covarianceEstAtPoint} estimates the covariance of a single vector $v$ by ratioing the set covariance of $\Xi$ to the set covariance of of observation window. Set covariance is calculated by intersecting a set with a translated copy of itself.

The others \code{covarianceMapEst_direct} and \code{covarianceRACS} estimate the covariance function at multiple $v$, returning 2D arrays of estimates.

\code{covarianceMapEst_direct} estimates covariance on a regular grid using the resolution of $\Xi_{obs}$. The regular grid extends to maxXshiftdistance and maxYshiftdistance in the x and y componenet respectively. It uses covarianceEstAtPoint to estimate the covariance at each of these points.

\code{covarianceRACS} is much faster (500 times faster in one comparison) because it uses Fourier Transforms. It uses SpatStat's setcov() function (which contains the Fourier Transform trick) to estimate set covariance over a window.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{comp1 }{An estimate (assuming stationarity of $\Xi$) that two points separated by $v$ will be in $\Xi$.}
\item{comp2 }{Denominator - The set covariance of the boundary}
\item{comp3 }{Numerator - The set covariance of $\Xi_{obs}}

For \code{covarianceEstAtPoint} these are single numerical values; for \code{covarianceMapEst_direct} they are matrices; for \code{covarianceRACS} they are objects of SpatStat's \code{im} class.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kassel Hingee
%%  ~~who you are~~
}
\note{
The name of this function is probably a bit confusing. Perhaps call it two-point coverage probability, and the other functions `two-point coverage probability' functions.

Double check that I'm calculating set covariance directly properly (I don't think there is a need to reflect the set, but maybe I got things wrong)

I should use "window" instead of "boundary"
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library("spatstat")

data(heather)
XiOWIN <- heather$coarse
windowOWIN <- owin(xrange=heather$coarse$xrange,yrange=heather$coarse$yrange)

coverageProb <- covpest(XiOWIN,windowOWIN)

covariancePt <- covarianceEstAtPoint(XiOWIN,windowOWIN,c(3,5))

covarianceDirectEst <- covarianceMapEst_direct(XiOWIN,windowOWIN,1,1)
filled.contour(covarianceDirectEst$covariance)


covarianceFcn <- covarianceRACS(XiOWIN,windowOWIN)
plot(covarianceFcn$covariance)
plot(covarianceFcn$covariance - coverageProb*coverageProb)


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
