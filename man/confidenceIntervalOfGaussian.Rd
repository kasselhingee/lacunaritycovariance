\name{confidenceIntervalOfGaussian}
\alias{confidenceIntervalOfGaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals from confidence levels.
%%  ~~function to do ... ~~
}
\description{
Gives a symmetric confidence interval for an estimator with Gaussain distribution (with mean 0).
This function probably already exists.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
confidenceIntervalOfGaussian(confidenceLevel, sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{confidenceLevel}{
  The desired probability (confidence) that the true value is contained in the returned range. 
   It should be a value between 0 and 1.
%%     ~~Describe \code{confidenceLevel} here~~
}
  \item{sd}{
  The standard deviation of the estimator.
%%     ~~Describe \code{sd} here~~
}
}
\details{
Just a wrapper of the quantile function qnorm()
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector (length = 2) of numerical values. First entry is the minimum of the range, second entry the maximum of the range.
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
 This function's behaviour is probably contained within a standard R function, but I figured this is good practise for me anyway.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
confidenceIntervalOfGaussian(0.95,1)
}