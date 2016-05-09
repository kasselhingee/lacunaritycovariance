## tests for Boolean simulations of deterministic disc code


library(stationaryracsinference, quietly = TRUE)

#check that its producing the same pattern
#Boolean model with discs of radius 10.
#The intensity has been chosen such that the true coverage fraction is very close to 0.5.
discr <- 10
w <- owin(xrange=c(0,100),c(0,100))
lambda <- 2.2064E-3
xi <- simulateBooleanDetermDiscs(lambda,discr,w,seed=6549)
coveragefrac(xi,w)

#theoretical coveragefrac
truecoveragefrac <- booldetermdiscs_truecoveragefrac(lambda,discr)

## thcovarDeterministicDiscs check
thcovariance <- thcovarDeterministicDiscs(xrange=c(-100,100),yrange=c(-100,100),eps=c(1,1),0.0022064,10)

## thspecdensAtOrigin check
thspecdens_origin <- thspecdensAtOrigin(lambda,discr)
thspecdens_origin

## use covariance to form an (approximate) confidence interval
setcovB <- setcov(w)
integrand <- eval.im((thcovariance-truecoveragefrac^2)*setcovB,harmonize = TRUE)
exactvariance <- (1/(area.owin(w))^2)*sum(integrand)*integrand$xstep*integrand$ystep
##assume that window large enough that distribution of estimator is Gaussian
##the exact variance is as calculated. Below is a 0.05 conf interval
confint0.05 <- c(qnorm(0.025, mean = truecoveragefrac, sd=sqrt(exactvariance)),
                 qnorm(1-0.025, mean = truecoveragefrac, sd=sqrt(exactvariance)))

##check that a randomly simulated value is in this range
xi <- simulateBooleanDetermDiscs(lambda,discr,w)
(coveragefrac(xi,w) < confint0.05[2]) & (coveragefrac(xi,w) > confint0.05[1])

