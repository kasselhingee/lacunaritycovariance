\name{exactvarP}
\alias{exactvarP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exact variance of coverage probability -- under construction
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
exactvarP(XiOWIN, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XiOWIN}{
%%     ~~Describe \code{XiOWIN} here~~
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
function (XiOWIN, w) 
{
    Xiinside <- intersect.owin(XiOWIN, w)
    p <- area.owin(Xiinside)/area.owin(w)
    setcovXi <- setcov(Xiinside)
    setcovB <- setcov(w)
    integrand <- eval.im(setcovXi - p^2 * setcovB, harmonize = TRUE)
    edgeValues = c(integrand[1, -1], integrand[-1, ncol(integrand)], 
        integrand[nrow(integrand), -ncol(integrand)], integrand[-nrow(integrand), 
            1])
    if (max(integrand, na.rm = TRUE) > diff(range(integrand, 
        na.rm = TRUE)) * 1e-05) {
        warning("covariance isn't uniformly close to p^2 at boundary\n")
        cat("max difference between C(v) and p^2 on boundary is ", 
            max(integrand, na.rm = TRUE), sep = "")
    }
    return((1/(area.owin(w))^2) * sum(integrand) * integrand$xstep * 
        integrand$ystep)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
