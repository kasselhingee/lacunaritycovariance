#functions for simulating general germ grain RACS

#randomly places grains with from a library of grains
#grainlib must a solist, pp a point pattern
placegrainsfromlib <- function(pp,grainlib,replace=TRUE,prob=NULL){
  grains <- sample(grainlib,size=pp$n,replace=replace,prob=prob) 
  pointlocations <- cbind(X=pp$x,Y=pp$y)
  pointlocations <- split(cbind(pointlocations),row(pointlocations)) #split matrix into a list of the rows
  shiftedgrains <- as.solist(mapply(shift.owin,grains,vec=pointlocations,SIMPLIFY=FALSE))
  placedgrains<- union.owin(shiftedgrains)  
  return(placedgrains)
}