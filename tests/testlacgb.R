suppressPackageStartupMessages(library(racssummfuncs))

img <- as.im(heather$coarse,eps=heather$coarse$xstep, na.replace=0)
sidel <- c(2.2)
lac <- lacgb(img,sidel,inclraw = TRUE)
lac$MVL
lac$raw
