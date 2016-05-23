## tests for Boolean simulations of deterministic disc code


library(stationaryracsinference, quietly = TRUE)

##test that seed option is keeping things reproducible
w <- owin(xrange=c(0,10),yrange=c(0,10))
xi <- rboollognormdiscs(w,2,1,-1,0.5,seed=36)
area.owin(xi)

