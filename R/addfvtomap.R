#' @title addfvtomap
#' @export addfvtomap
#'
#' @description A function for plotting fv objects onto large spatial map figures. This allows localised summary functions to be plotted on the spatial region that they describe.
#' 
#' @param fvobj Is the fv function to plot
#' @param mapxlim is a vector of length 2 giving the x-coordinates of the region to use to plot \code{fvobj} in the units of the current figure.
#' @param mapylim Same as mapxlim except for the y-coordinates.

#' @return Nothing. Adds a plot of the default function value in \code{fvobj} to the current figure within the region given by mapxlim and mapylim.

#' @example


addfvtomap <- function(fvobj, mapxlim, mapylim){
  fvobj.ymin <- range.fv(fvobj)[[1]]
  fvobj.ymax <- range.fv(fvobj)[[2]]
  fvobj.xmin <- attributes(fvobj)$alim[[1]]
  fvobj.xmax <- attributes(fvobj)$alim[[2]]
  p.ymin <- mapylim[[1]]
  p.ymax <- mapylim[[2]]
  p.xmin <- mapxlim[[1]]
  p.xmax <- mapxlim[[2]]
  yscale <- (p.ymax - p.ymin) / (fvobj.ymax - fvobj.ymin)
  xscale <- (p.xmax - p.xmin)/(fvobj.xmax - fvobj.xmin)
  plot.fv(add = TRUE, fvobj,
          p.ymin + yscale * (.y - fvobj.ymin) ~ p.xmin + xscale * (.x - fvobj.xmin),
          xlim = mapxlim,
          ylim = mapylim,
          col = "green"
          )
  arrows(mapxlim[[1]], mapylim[[1]], mapxlim[[2]], mapylim[[1]], lend = 1, angle = 90, length = 0)
  arrows(mapxlim[[1]], mapylim[[1]], mapxlim[[1]], mapylim[[2]], lend = 1, angle = 90, length = 0)
  text(mean(mapxlim), y = mapylim[[1]], pos = 1, labels = format(fvobj.xmax - fvobj.xmin, digits = 2))
  text(mapxlim[[1]], y = mean(mapylim), pos = 2, labels = format(fvobj.ymax - fvobj.ymin, digits = 2))
}
