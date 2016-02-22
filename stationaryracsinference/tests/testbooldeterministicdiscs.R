## tests for Boolean simulations of deterministic disc code


library(stationaryracsinference, quietly = TRUE)

## thcovarDeterministicDiscs check
thcovarDeterministicDiscs(xrange=c(-10,10),yrange=c(-3,3),eps=c(1,1),0.0022064,10)

## thspecdensAtOrigin check
thspecdens_origin <- thspecdensAtOrigin(0.0022064,10)
thspecdens_origin

## quasithspecdens - check at origin
specdens <- quasithspecdens(0.0022064,10)
specdens[round(dim(specdens)[2]/2),round(dim(specdens)[1]/2)]

