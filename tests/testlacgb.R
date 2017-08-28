suppressPackageStartupMessages(library(racssummfuncs))

img <- as.im(heather$coarse,eps=heather$coarse$xstep, na.replace=0)
sidel <- c(2.2)
lac.wraw <- lacgb(img,sidel,inclraw = TRUE)
lac.wraw$MVL
lac.wraw$raw

lac <- lacgb(img,sidel)
lac$MVL
