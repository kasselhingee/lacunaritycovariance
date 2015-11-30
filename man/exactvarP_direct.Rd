\name{exactvarP_direct}
\alias{exactvarP_direct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exact variance of coverage probability ... under construction
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
exactvarP_direct(pest, covarianceMap, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pest}{
%%     ~~Describe \code{pest} here~~
}
  \item{covarianceMap}{
%%     ~~Describe \code{covarianceMap} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pest, covarianceMap, w) 
{
    edgeValues = c(covarianceMap$covariance[1, -1], covarianceMap$covariance[-1, 
        ncol(covarianceMap$covariance)], covarianceMap$covariance[nrow(covarianceMap$covariance), 
        -ncol(covarianceMap$covariance)], covarianceMap$covariance[-nrow(covarianceMap$covariance), 
        1])
    if (max(abs(edgeValues - pest^2), na.rm = TRUE) > diff(range(covarianceMap$covariance, 
        na.rm = TRUE)) * 1e-05) {
        warning("covariance map doesn't reached p^2 at boundary")
    }
    return(1/(area.owin(w)^2) * vectorSteps * sum(covarianceMap$numerator - 
        covarianceMap$denominator * pest^2, na.rm = TRUE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
