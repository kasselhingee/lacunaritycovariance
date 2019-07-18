#' @title Test if an \code{im} object is a binary map
#' @export isbinarymap
#' @import spatstat
#' @param xi an image object
#' @param requiretrue Logical. If TRUE then \code{isbinarymap} will error if xi is NOT a binary map.
#' @description
#' Tests whether \code{xi} is a binary map. 
#' The pixel values must be of logical type (\code{TRUE}, \code{FALSE} and \code{NA} only), or numerical (1, 0 or \code{NA}).

#' @return Logical value. \code{TRUE} if \code{xi} is a binary map. Otherwise \code{FALSE}.
#' If \code{requiretrue = TRUE} and \code{xi} is not a binary map then an error will occur.

#' @examples
#' # The following return TRUE
#' isbinarymap(as.im(heather$coarse, na.value = 0))
#' isbinarymap(as.im(heather$coarse, na.value = FALSE, value = TRUE))
#' 
#' # The following return FALSE
#' isbinarymap(as.im(heather$coarse, na.value = 0.2, value = 1))
#' isbinarymap(as.im(heather$coarse, na.value = 0, value = 1.5))

isbinarymap <- function(xi, requiretrue = FALSE) {
    uvals <- unique(as.list(as.matrix(xi)))
    if ( !all(  (uvals %in% c(0, 1)) | is.na(uvals))  && 
             !all((uvals %in% c(FALSE, TRUE, NA)) | is.na(uvals)) ) {
        isbinary <- FALSE
    } else {
        isbinary <- TRUE
    }
    if (requiretrue && !isbinary) {stop("Input xi has values other than 0, 1 or NA")}
    return(isbinary)
}
