suppressPackageStartupMessages(library(stationaryracsinference))

img <- as.im(heather$coarse,eps=heather$coarse$xstep)
bandwidths <- c(1.1)
lac <- lacgb(img,bandwidths)
lac$RS
lac$nobord
