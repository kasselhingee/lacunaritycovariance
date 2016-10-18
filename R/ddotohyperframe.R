#' @title trelliscope ddo to spatstat hyperframe conversion
#' @export ddotohyperframe
#' 
#' @description Extracts data from tessera's ddo objects and puts them into a hyperframe.
#' 
#' @param inddo Input ddo object
#' @param slots A list of names of slots to extract into a hyperframe - each slot becomes a column in the hyperframe. 
#' If not supplied then function will try to extract all slots using the first key-value pair as an example.
#' @return A hyperframe, each column corresponding to a slot.
#' 
###########################################################
#' 
#' @examples
#' applyestimators <- function(obs,w){
  #' #coverage fraction
  #' coveragefrac <- coveragefrac(obs,w)
  
  #' #covariance
  #' covariances <- covariance(obs,w)
  #' 
  #' #Sph Contact Functions. 
  #' sphContact <- list(xiHest = Hest(obs),
  #'                    notxiHest = Hest(complement.owin(obs))) 
  #' 
  #' #contagion estimates
  #' sphContactHarm <-  harmonise(xiHest = sphContact[[1]],notxiHest = sphContact[[2]])
  #' sphcontactcontag <- list(
  #'   contagionEst = contagSphCont(sphContactHarm$xiHest$km,
  #'                                sphContactHarm$notxiHest$km,
  #'                                coveragefrac),
  #'   r = sphContactHarm$xiHest$r)
  #' 
  #' return(list(
  #'   coveragefrac = coveragefrac,
  #'   covariance = covariances,
  #'   sphcontactxi = sphContact$xiHest,
  #'   sphcontactnotxi = sphContact$notxiHest,
  #'   sphcontactcontag = sphcontactcontag
  #' ))
#' }
#' 
#' estimates <- datadr::kvPairs(
#'      datadr::kvPair("medium",applyestimators(heather$medium,Frame(heather$medium))),
#'      datadr::kvPair("coarse",applyestimators(heather$coarse,Frame(heather$coarse)))
#'      )
#' estimatesddo <- datadr::ddo(estimates)
#' 
#' hframe <- ddotohyperframe(estimatesddo,slots = c("coveragefrac","sphcontactxi"))
#' hframe <- ddotohyperframe(estimatesddo)
#' plot(hframe[,"sphcontactxi",drop=TRUE], main = "Spherical Contact Distribution:\n Diggle's Heather Data")
#' 
ddotohyperframe <- function(inddo,slots="all"){
  if (length(slots)==1 && slots == "all"){
    warning("`all' option uses names of the first kv pair, if the values of each pair are different
            then this wont capture everyting. A more robust method is to specify the slots you are interested in.")
    slots <- names(inddo[[1]][[2]]) #names of object stored in the value of the first kvpair
  } 
  hframe <- hyperframe()
  for (slotn in slots){
    kvpairslist <- datadr::drLapply(inddo, function(v) v[[slotn]],
                                    combine=datadr::combCollect,
                                    params=c(slotn = slotn))  #returns a list of kvpairs
    #note that the params field doesn't seem to be working
    # - changing the parameter name to slotm makes the function error
    
    valueslist <- lapply(kvpairslist,"[[","value")
    names(valueslist) <- lapply(kvpairslist,"[[","key")
    
    #don't do anything if all of length one then unlist! Objects like windows are all one length to. The trick will be to save them in the ddo as non-nested objects
    hframe <- cbind(hframe,valueslist)
  }
  row.names(hframe) <- names(valueslist)
  names(hframe) <- slots
  return(hframe)
}

#' addTransform(inddo, function(v) v$lpicontag[[1]])
#' 
