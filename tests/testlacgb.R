suppressPackageStartupMessages(library(stationaryracsinference))

img <- as.im(heather$fine,eps=heather$fine$xstep)
bandwidths <- 2
lac <- lacgb(img,bandwidths)
abs(lac$RS-0.07396426) < 1E-6
abs(lac$nobord-0.279406) < 1E-6
