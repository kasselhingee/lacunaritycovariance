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
#' @return An fv object containing the pointwise mean, variance and maxima and minima.

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
  nacount <- eval.fv(0 * is.na(A), envir = list(A = object[[1]]), dotonly = FALSE, relabel = FALSE)
  if (na.rm){
    #count number of NA values
    nas <- Map(function(a) eval.fv(is.na(a), dotonly = FALSE, relabel = FALSE), object)
    nacount <- Reduce(Add.fv, nas)
    #remove NA vals from object
    object <- Map(function(a) replace(a, is.na(as.data.frame(a)), 0), object)
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
                        lablprefix="")
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
Add.fv <- function(A,B){ force(A); force(B); eval.fv(A+B, dotonly = FALSE, relabel=FALSE) }

Pmax.fv <- function(A, B){force(A); force(B); eval.fv(pmax(A, B), dotonly = FALSE, relabel = FALSE)}
Pmin.fv <- function(A, B){force(A); force(B); eval.fv(pmin(A, B), dotonly = FALSE, relabel = FALSE)}


#this is a draft function - it doesn't operate for every function value in an fv object!
summary_fvlist_fivenum <- function(object, ...){
  object <- harmonise.fv(object)
  #converting into data.frames
  ynames <- fvnames(object[[1]], a = ".")
  xname <- fvnames(object[[1]], a = ".x")
  yvals <- lapply(object, function(x) as.matrix(x[, ynames[[1]], drop = TRUE]))
  yvals.m <- do.call(cbind, args = yvals) #each column is an fv object, each row is a unique argument value
  ptwisesum <- apply(yvals.m, MARGIN = 1, FUN = fivenum) #each column is a unique argument value, each row is a quartile thing
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
