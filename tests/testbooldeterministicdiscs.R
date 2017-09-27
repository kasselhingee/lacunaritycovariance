## tests for Boolean simulations of deterministic disc code


suppressPackageStartupMessages(library(stationaryracsinference))

#check that its producing the same pattern
#Boolean model with discs of radius 10.
#The intensity has been chosen such that the true coverage fraction is very close to 0.5.
discr <- 10
w <- owin(xrange=c(0,100),c(0,100))
lambda <- 2.2064E-3
xi <- rbdd(lambda,discr,w,seed=6549)
coveragefrac(xi,w)

#theoretical coveragefrac
truecoveragefrac <- bdd.coverageprob(lambda,discr)

## bdd.covar check
thcovariance <- bdd.covar(xrange=c(-100,100),yrange=c(-100,100),eps=c(1,1),0.0022064,10)

## bdd.specdensAtOrigin check
thspecdens_origin <- bdd.specdensAtOrigin(lambda,discr)
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
xi <- rbdd(lambda,discr,w)
(coveragefrac(xi,w) < confint0.05[2]) & (coveragefrac(xi,w) > confint0.05[1])

