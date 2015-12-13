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
HestInPolygonalWindow(X, boundaryIn = owin(X$xrange, X$yrange),\
r = NULL, breaks = NULL, ..., correction = c("km", "rs", "han"),\
conditional = TRUE)
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

XiOWIN <- heather$coarse
windowOWIN <- Frame(heather$coarse)

Z <- HestInPolygonalWindow(XiOWIN,windowOWIN)
plot(Z)
plot(add=TRUE,Hest(XiOWIN))#since window is rectangular these should be the same

}
