suppressPackageStartupMessages(library(racssummfuncs))

covar <- covariance(heather$coarse)
p <- area(heather$coarse)/area(Frame(heather$coarse))
sidelengths <- 2.2
lac <- lac(sidelengths,covar,p)
lac$MVL
