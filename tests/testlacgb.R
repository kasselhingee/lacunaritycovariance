suppressPackageStartupMessages(library(stationaryracsinference))

img <- as.im(heather$coarse,eps=heather$coarse$xstep)
sidel <- c(2.2)
lac <- lacgb(img,sidel)
lac$RS
lac$raw
