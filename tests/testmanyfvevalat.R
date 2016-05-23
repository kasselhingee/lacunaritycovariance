
library(stationaryracsinference, quietly = TRUE)

fv1 <- Hest(heather$coarse)
fv2 <- Hest(complement.owin(heather$coarse))
fvlist <- list(fv1,fv2) 

fvfunc01 <- as.function.fv(fv1, value="km")
manyfvevalat(fvlist,0.19,value="km")[[1]]==fvfunc01(0.19)

manyfvevalat(fvlist,0.19,value="km")

