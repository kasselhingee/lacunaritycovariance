#' @title Estimates the integral of image values within a radius r.
#' @export integratebyradius 
#' 
#' @description Estimates the integral of image values within a radius r from the centre.
#' @return  An fv object
#' @param centre The centre location
#' @param values A pixel image to be integrated

#' @examples 
#' twptprob <- covariance(heather$coarse,Frame(heather$coarse))
#' kfcn <- integratebyradius(c(0,0),twptprob)
#' plot(kfcn)
#' 
#' 
integratebyradius <- function(centre, values){
  #this will be a lot faster if I do the whole thing in integer arithmetic, 
  #which is definitely possible for images with equal x and y steps
  sqdist <- outer((values$yrow-centre[[2]])^2, #different yrow here corresponds to different rows in sqdist
        (values$xcol-centre[[1]])^2, #different xcol here correpond to different cols in sqdist
        FUN ="+")
  #each col in above sqdist corresponds to the same xcol
  #each row in above sqdist corresponds to the same yrow
  valuesm <- as.matrix(values) #a row in valuesm correspond to a single y value
  
  #turn into long lists and sort
  sqdist <- as.vector(sqdist)
  valuesm <- as.vector(valuesm)
  perm <- order(sqdist,decreasing=FALSE)
  sqdist <- sqdist[perm]
  valuesm <- valuesm[perm]
  
  #cut out everything above an NA value
  valuesm <- valuesm[1:which.min(!is.na(valuesm))]
  sqdist <- sqdist[1:length(valuesm)]
  
  #compute shortest distance to boundary from centre - stop summing after that happens.
  distfunc <- distfun.owin(Frame(values), discretise = FALSE, invert=TRUE)
  shdist <- distfunc(centre[[1]],centre[[2]])
  sqshdist <- shdist*shdist
  csumvals <- cumsum(valuesm[sqdist<sqshdist])
  sqdist <- sqdist[sqdist<sqshdist]
  
  #now to only take last sqdist value in lists
  keepind <- duplicated(sqdist,fromLast = TRUE) #indox of last of each unique sqdist
  rvals <- sqrt(sqdist[keepind])
  intwithinrdist <- csumvals[keepind]*values$xstep*values$ystep
  cumint <- fv(data.frame(r = rvals,
                          y = intwithinrdist),
               valu = "y",
               desc=c("radius",
                      "Estimated integral of values within radius r"),
               alim = range(rvals)
               )
  return(cumint)
}
