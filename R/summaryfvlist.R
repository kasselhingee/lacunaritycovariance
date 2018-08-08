#' @title Summarise a list fv objects
#' @export summary.fvlist
#' @method summary fvlist
#' @description
#' This function assumes that fv objects are each realisations of the same (stochastic) object. 
#' It returns pointwise summaries such as observed sample mean and sample variance.

#' @author Kassel Hingee

#' @param  object A list of fv objects
#' @param ...  Ignored.
#' @param na.rm If TRUE NA values in the fv object will be ignored 
#' (i.e. NA values will be removed from all summations and the population size will decrease by 1 for each NA value).
#' @return An fv object containing the pointwise mean, pointwise sample variance, pointwise maximum and pointwise minimum.

#' @examples
#' obspatterns <- rpoispp(10, nsim = 10)
#' object <- lapply(obspatterns, Hest, W = Frame(obspatterns[[1]]))
#' summ <- summary.fvlist(object)
#' 
#' #test with NA vals
#' object[[1]]$km[1:46] <- NA
#' summ <- summary.fvlist(object, na.rm = TRUE)

summary.fvlist <- function(object, ..., na.rm = FALSE){
  object <- harmonise.fv(object)
  n <- length(object)
  nacount <- with.fv(object[[1]], 0 * ., fun = TRUE)
  if (na.rm){
    #count number of NA values
    nas <- Map(function(a) eval.fv(is.na(a), dotonly = FALSE, relabel = FALSE), object)
    nacount <- Reduce(Add.fv, nas)
    
    #use narm versions of Add.fv, Pmax.fv and Pmin.fv functions
    Add.fv <- Add.fv.narm
    Pmax.fv <- Pmax.fv.narm
    Pmin.fv <- Pmin.fv.narm
  } else {#set nakeep versions of  Add.fv, Pmax.fv and Pmin.fv etc
    Add.fv <- Add.fv.nakeep
    Pmax.fv <- Pmax.fv.nakeep
    Pmin.fv <- Pmin.fv.nakeep
  }
  
  
  #compute mean
  meanY <- Reduce(Add.fv, object) #compute sum
  meanY <- eval.fv(meanY/(n - nacount), dotonly = FALSE, relabel = FALSE) #divide by number of non-NA vals
  
  #compute variance
  sumY2 <- Reduce(Add.fv, lapply(object, Square.fv)) #use NA removed version
  varY <- eval.fv( (sumY2 - (n - nacount) *( meanY^2))/(n - nacount - 1), dotonly = FALSE, relabel = FALSE)
  varY <- eval.fv(pmax.int(0, varY), dotonly = FALSE, relabel = FALSE)
  maxY <- Reduce(Pmax.fv, object)
  minY <- Reduce(Pmin.fv, object)
  ## tweak labels of main estimate
  attributes(meanY) <- attributes(varY) <- attributes(maxY) <- attributes(minY) <- attributes(vanilla.fv(object[[1]]))
  attributes(minY) <- attributes(vanilla.fv(object[[1]]))
  meanY <- prefixfv(meanY,
                        tagprefix="mean",
                        descprefix="mean ",
                        lablprefix="bold(mean)~")
    ## tweak labels of variance terms
  varY <- prefixfv(varY,
                         tagprefix="var",
                         descprefix="sample variance of ",
                         lablprefix="bold(var)~")
  maxY <- prefixfv(maxY,
                         tagprefix="max",
                         descprefix="maximum ",
                         lablprefix="bold(max)~")
  minY <- prefixfv(minY,
                         tagprefix="min",
                         descprefix="minimum ",
                         lablprefix="bold(min)~")
    ## glue together
    result <- Reduce(bind.fv, list(meanY, varY, maxY, minY))
    # set default plotting lines
    fvnames(result, ".") <- c(fvnames(meanY, "."), fvnames(maxY, "."), fvnames(minY, "."))
    attr(result, "fmla") <- ". ~ .x"
  return(result)
}

#Handy functions copied from spatstat
Square.fv <- function(A) { force(A); eval.fv(A^2, dotonly = FALSE, relabel=FALSE) }
Add.fv.nakeep <- function(A,B){ force(A); force(B); eval.fv(A+B, dotonly = FALSE, relabel=FALSE) }
Add.fv.narm <- function(A,B){
  force(A)
  force(B)
  #do addition that keeps NAs
  out <- Add.fv.nakeep(A, B)
  
  #replace any NAs with values from A or B if they exist
  fillfromA <- is.na(as.data.frame(out)) & (!is.na(as.data.frame(A)))
  out <- replace(out, fillfromA, as.data.frame(A)[fillfromA])
  fillfromB <- is.na(as.data.frame(out)) & (!is.na(as.data.frame(B)))
  out <- replace(out, fillfromB, as.data.frame(B)[fillfromB])
  if (any(fillfromA & fillfromB)) {stop("Add.fv.narm unable to remove the some NA values.")}
  return(out)
  }

Pmax.fv.nakeep <- function(A, B){force(A); force(B); eval.fv(pmax(A, B, na.rm = FALSE), dotonly = FALSE, relabel = FALSE)}
Pmax.fv.narm <- function(A, B){force(A); force(B); eval.fv(pmax(A, B, na.rm = TRUE), dotonly = FALSE, relabel = FALSE)}
Pmin.fv.nakeep <- function(A, B){force(A); force(B); eval.fv(pmin(A, B, na.rm = FALSE), dotonly = FALSE, relabel = FALSE)}
Pmin.fv.narm <- function(A, B){force(A); force(B); eval.fv(pmin(A, B, na.rm = TRUE), dotonly = FALSE, relabel = FALSE)}

#this is a draft function - it doesn't operate for every function value in an fv object!
summary_fvlist_fivenum <- function(object, ...){
  object <- harmonise.fv(object)
  #converting into data.frames
  ynames <- fvnames(object[[1]], a = ".")
  xname <- fvnames(object[[1]], a = ".x")
  yvals <- lapply(object, function(x) as.matrix(x[, ynames[[1]], drop = TRUE]))
  yvals.m <- do.call(cbind, args = yvals) #each column is an fv object, each row is a unique argument value
  ptwisesum <- apply(yvals.m, MARGIN = 1, FUN = stats::fivenum) #each column is a unique argument value, each row is a quartile thing
  fvdata <- as.data.frame(cbind(object[[1]][, xname, drop = TRUE], t(ptwisesum)))
  names(fvdata) <- c(xname, "min", "Q1", "median", "Q3", "max")
  fivenumfv <- fv(fvdata,
     argu = xname,
     ylab = attr(object[[1]], "ylab", exact = TRUE), 
     valu = "median",
     fmla = ". ~ r",
     alim = c(min(fvdata[, xname]), max(fvdata[, xname])),
     labl = names(fvdata),
     desc = c(attr(object[[1]], "desc", exact = TRUE)[[1]], 
              "minimum",
              "1st quartile",
              "median",
              "3rd quartile",
              "maximum"),
     unitname = unitname(object[[1]]),
     fname = paste("Five Number Summary of", attr(object[[1]], "fname"))
  )
  return(fivenumfv)
}
