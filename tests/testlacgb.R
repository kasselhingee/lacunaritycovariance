library(stationaryracsinference, quietly = TRUE)

data(balcattapark_coarse)
img <- as.im(balcattapark_coarse$vegmask)
bandwidths <- 5*0.8
lac <- lacgb(img,bandwidths)
abs(lac$RS-1.701478) < 1E-6
abs(lac$nobord-1.776974) < 1E-6
