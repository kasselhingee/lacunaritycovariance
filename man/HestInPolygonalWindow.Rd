\name{HestInPolygonalWindow}
\alias{HestInPolygonalWindow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate the empty space function from a non-rectangular observation window
%%  ~~function to do ... ~~
}
\description{A hack of SpatStat's Empty Space estimator.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HestInPolygonalWindow(X, boundaryIn = owin(X$xrange, X$yrange), r = NULL, breaks = NULL, ..., correction = c("km", "rs", "han"), conditional = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{boundaryIn}{
%%     ~~Describe \code{boundaryIn} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{correction}{
%%     ~~Describe \code{correction} here~~
}
  \item{conditional}{
%%     ~~Describe \code{conditional} here~~
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
Mostly others...
%%  ~~who you are~~
}
\note{The difference to SpatStat Hest(), is that this accepts an observation window argument and it calculates the censoring distance using this window. 
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(spatstat)

data(heather)
XiOWIN <- heather$coarse
windowOWIN <- owin(xrange=heather$coarse$xrange,yrange=heather$coarse$yrange)

Z <- HestInPolygonalWindow(XiOWIN,windowOWIN)
plot(Z)
plot(add=TRUE,Hest(XiOWIN))#since window is rectangular these should be the same

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, boundaryIn = owin(X$xrange, X$yrange), r = NULL, 
    breaks = NULL, ..., correction = c("km", "rs", "han"), conditional = TRUE) 
{
    if (!(is.ppp(X) || is.psp(X) || is.owin(X))) 
        stop("X should be an object of class ppp, psp or owin")
    if (is.null(correction)) 
        correction <- c("rs", "km", "cs")
    correction <- pickoption("correction", correction, c(none = "none", 
        raw = "none", border = "rs", rs = "rs", KM = "km", km = "km", 
        Kaplan = "km", han = "han", Hanisch = "han", best = "km"), 
        multi = TRUE)
    corxtable <- c("km", "rs", "han", "none")
    corx <- as.list(corxtable \%in\% correction)
    names(corx) <- corxtable
    X <- intersect.owin(X, boundaryIn)
    D <- distmap(X)
    B <- distmap.owin(boundaryIn, dimyx = X$dim, invert = TRUE)
    W <- as.owin(D)
    dmax <- summary(D)$max
    breaks <- handle.r.b.args(r, breaks, W, NULL, rmaxdefault = dmax)
    rval <- breaks$r
    dist <- as.vector(as.matrix(D))
    bdry <- as.vector(as.matrix(B))
    ok <- !is.na(dist) && !is.na(bdry)
    dist <- dist[ok]
    bdry <- bdry[ok]
    if (is.owin(X)) {
        pos <- (dist > 0)
        inWindow <- (bdry >= 0)
        areafraction <- 1 - (sum(pos & inWindow)/sum(inWindow))
        dist <- dist[pos & inWindow]
        bdry <- bdry[pos & inWindow]
    }
    d <- (dist <= bdry)
    o <- pmin.int(dist, bdry)
    Z <- censtimeCDFest(o, bdry, d, breaks, KM = corx$km, RS = corx$rs, 
        HAN = corx$han, RAW = corx$none, han.denom = if (corx$han) 
            eroded.areas(W, rval)
        else NULL, tt = dist)
    if (is.owin(X)) {
        if (conditional) {
            if (corx$km) 
                Z$km <- condition(Z$km)
            if (corx$rs) 
                Z$rs <- condition(Z$rs)
            if (corx$han) 
                Z$han <- condition(Z$han)
            if (corx$none) 
                Z$raw <- condition(Z$raw)
        }
        else {
            if (corx$km) 
                Z$km <- reconstitute(Z$km, areafraction)
            if (corx$rs) 
                Z$rs <- reconstitute(Z$rs, areafraction)
            if (corx$han) 
                Z$han <- reconstitute(Z$han, areafraction)
            if (corx$none) 
                Z$raw <- reconstitute(Z$raw, areafraction)
        }
    }
    Z <- rebadge.fv(Z, substitute(H(r), NULL), "H")
    unitname(Z) <- unitname(X)
    return(Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
