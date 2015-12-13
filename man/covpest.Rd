\name{covpest}
\alias{covpest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate the coverage probability of a stationary RACS \eqn{P(o)\in \Xi}.
%%  ~~function to do ... ~~
}
\description{Estimate is simply the fraction of \eqn{\Xi} in the observed window
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
covpest(Xi, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xi}{An observation of a RACS.
%%     ~~Describe \code{Xi} here~~
}
  \item{w}{The window of observation (not necessarily rectangular)
%%     ~~Describe \code{w} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{An estimate of the coverage probability
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
XiOWIN <- heather$coarse
windowOWIN <- Frame(heather$coarse)

coverageProb <- covpest(XiOWIN,windowOWIN)

#vegetation map of balcatta park
data(balcattapark_coarse)
coverageProb <- covpest(balcattapark_veg_coarse,balcattapark_boundary)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
