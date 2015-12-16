#function for simulating a Boolean Model of discs with lognormal radius

rboollognormdiscs <- function(window,bufferdist,lambda,meanlog,sdlog){
  #have to simulate in a much larger area than the observation window (because grains with centres outside the window should still be observed)
  wsim <- Frame(dilation(window,bufferdist)) #i reckon faster to use rectangular region (the non-rectangular probably simulates in a rectangular region and then rejects anyway)
  pp <- rpoispp(lambda,win=wsim,nsim=1,drop=TRUE) #prepare a random radius for each point
 
  radius <- rlnorm(pp$n,meanlog=meanlog,sdlog=sdlog) #prepare a random radius for each point
   
  #calculating grains
  pointlocations <- cbind(X=pp$x,Y=pp$y)
  pointlocations <- split(cbind(pointlocations),row(pointlocations)) #split matrix into a list of the rows
  grains <- mapply(disc,radius = radius,centre=pointlocations,SIMPLIFY=FALSE) #calculate grains with their locations
  
  #take union of all grains
  xisim <- union.owin(as.solist(grains))
  
  xi <- intersect.owin(xisim,window)
  return(xi) 
}