suppressPackageStartupMessages(library(stationaryracsinference))

covar <- covariance(heather$coarse)
p <- area(heather$coarse)/area(Frame(heather$coarse))
bandwidths <- 2
lac <- lac(bandwidths,covar,p)
abs(lac-0.06496912) < 1E-6
