#' @title Plot multiple fv objects of similar ilk
#' @export manylines.fv
#' 
#' @description Plots multiple fv objects on the same window. 
#' @param fvlist A list of fv objects
#' @param xlim Plotting limits of the horizontal direction
#' @param ylim Plotting limits of the vertical direction
#' @param ... Plotting parameters to be passed to \code{title} and \code{lines}
#' @param xname the name of the value to associate with the x-axis
#' @param yname the name of the value to associate with the y-axis
#' @param add If true will add lines to current plot, otherwise will start new plot and will require xlim and ylim.
#' @return List of return values from lines (typically a list of NULL values) 
manylines.fv <- function(fvlist, xname, yname, ..., add=FALSE, xlim=NULL, ylim=NULL){
    if ((!add) & (is.null(xlim) | is.null(ylim))){stop("To create a new plot need xlim and ylim")}
    
    if (!add){
      plot.new()
      plot.window(xlim=xlim, ylim=ylim, ...)
      title(...)
      axis(1)
      axis(2)
    }
    
    if (is.null(names(fvlist))){names(fvlist) <- 1:length(fvlist)}
    
    xvals <- lapply(fvlist, "[[", xname) #extract lists of x values
    nullxvals <- as.logical(lapply(xvals,is.null))
    if (any(nullxvals)){
      warning(paste(yname, "is empty for ",paste(names(fvlist)[which(nullxvals)], collapse=", ")))
    }
    yvals <- lapply(fvlist, "[[", yname) #extract lists of y values
    nullyvals <- as.logical(lapply(yvals,is.null))
    if (any(nullyvals)){
      warning(paste(yname, "is empty for ",paste(names(fvlist)[which(nullyvals)], collapse=", ")))
    }
    out <- mapply(lines, xvals, yvals,
                  #col=rainbow(length(fvlist)),
                  ...)
    return(out)
}


#' @examples 
#' hestf <- Hest(heather$fine, W =Frame(heather$fine))
#' hestm <- Hest(heather$medium, W =Frame(heather$medium))
#' hestc <- Hest(heather$coarse, W =Frame(heather$coarse))
#' fvlist <- list(hestf,hestm,hestc)
#' names(fvlist) <- c("fine", "medium", "coarse")
#' manylines.fv(fvlist, xname="r", yname="km", xlim=c(0,0.35), ylim=c(0,1))
#'
#' ## Fancier plot
#' out <- manylines.fv(fvlist, xname="r", yname="km", xlim=c(0,0.35), ylim=c(0,1), 
#' lwd=c(1,2,3), main="Spherical Contact Distribution\n from different resolutions", col=rainbow(length(fvlist)))
