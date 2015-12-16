\name{rboollognormdiscs}
\alias{rboollognormdiscs}
\title{Simulates a Boolean model of discs with log normal disc radii and Poisson germs}
\description{
Simulates a Poisson point process and places discs of random radii around each point.
The radii are generated using a log normal distribution.
}
\usage{
rboollognormdiscs(window,bufferdist,lambda,meanlog,sdlog)
}
\arguments{
  \item{window}{An \code{owin} object specifying the desired simulation region}
  \item{bufferdist}{A distance to expand \code{window} so that discs with centres near \code{window} are also simulated}
  \item{lambda}{Intensity of the Poisson point process, passed to \code{\link{spatstat}{rpoispp}}.
   It could be either a single positive number, a function (x,y,...) or a pixel image 
   (they must be defined for the dilate window)}
  \item{meanlog}{for the distribution of radii. The logarithm of the distribution has mean equal to \code{meanlog}.}
  \item{sdlog}{For the distribution of radii. The logarithm of the distribution has standard deviation equal to \code{sdlog}}
}
\details{
The point process needs to be simulated in a larger region than the desired window to account for the possibility of discs that intersect the window, but have germs outside the window.

The point process is generated using spatstat's \code{\link[spatstat]{rpoispp}}.
}
\value{Returns an \code{owin} object cropped to \code{window}.}
\references{
[1] I. Molchanov, Statistics of the Boolean Model for Practitioners and Mathematicians. Wiley, 1997. (Chapter 8)

[2] **seminal paper - haven't read** Bindrich U. and Stoyan D. (1991) Stereology for pores in white bread: statistical analyses for the Boolean model by serial sections. J. Microscopy 162:231-239
}
\author{
Kassel Hingee
}
\section{Warning }{A good choice of bufferdist is required and it probably depends on the distribution of radii.} 

\examples{
w <- owin(xrange=c(0,10),yrange=c(0,10))
xi <- rboollognormdiscs(w,2,1,-1,0.5)

plot(w)
plot(xi,add=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
