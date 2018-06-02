#' @title Plot multiple fv objects of similar ilk
#' @export manylines.fv
#' @importFrom graphics plot.new plot.window title axis
#' 
#' @description Plots multiple fv objects on the same window. 
#' @param fvlist A list of fv objects
#' @param xlim Plotting limits of the horizontal direction
#' @param ylim Plotting limits of the vertical direction
#' @param ... Plotting parameters to be passed to \code{plot.window}, \code{title} and \code{\link{plot.fv}}
#' @param fmla an R language formula determining which variables or expressions are plotted.
#'  Either a formula object, or a string that can be parsed as a formula.
#'  It is passed to \code{plot.fv} - see \code{\link{plot.fv}} for details.
#' @param add If true will add lines to current plot, otherwise will start new plot and will require xlim and ylim.
#' @return List of return values from plot.fv 

#' @examples 
#' hestf <- Hest(heather$fine, W =Frame(heather$fine))
#' hestm <- Hest(heather$medium, W =Frame(heather$medium))
#' hestc <- Hest(heather$coarse, W =Frame(heather$coarse))
#' fvlist <- list(hestf,hestm,hestc)
#' names(fvlist) <- c("fine", "medium", "coarse")
#' # manylines.fv(fvlist, km~r)
#'
#' ## Fancier plot
#' # out <- manylines.fv(fvlist,  
#' # main="Spherical Contact Distribution\n from different resolutions", col=rainbow(length(fvlist)))
manylines.fv <- function(fvlist, fmla=NULL, ..., add=FALSE, xlim=NULL, ylim=NULL){
  if (!add && (is.null(xlim) || is.null(ylim))){#if not add and xlim ylim are null then
      limits <- mapply(plot.fv,fvlist,MoreArgs=list(fmla = fmla, limitsonly=TRUE), SIMPLIFY=FALSE)
      #get largest x limits
      if (is.null(xlim)){
        xlims <- lapply(limits,"[[","xlim")
        xmin <- min(unlist(lapply(xlims,"[[",1)))
        xmax <- max(unlist(lapply(xlims,"[[",2)))
        xlim <- c(xmin,xmax)       
      }

      #get largest y limits
      if (is.null(ylim)){
        ylims <- lapply(limits,"[[","ylim")
        ymin <- min(unlist(lapply(ylims,"[[",1)))
        ymax <- max(unlist(lapply(ylims,"[[",2)))
        ylim <- c(ymin,ymax)
      }
  }
  if (!add){ #if not adding then create new plot window
      plot.new()
      plot.window(xlim=xlim, ylim=ylim, ...)
      title(...)
      axis(1)
      axis(2)
  }
  
  if (is.null(names(fvlist))){names(fvlist) <- 1:length(fvlist)} #give fvlist names, for autolegend if the ... argument asks for it. 
  
  out <- mapply(plot.fv, fvlist, list(fmla=fmla), add=TRUE, ..., MoreArgs=list(xlim=xlim,ylim=ylim), SIMPLIFY = FALSE)
  return(out)
}
