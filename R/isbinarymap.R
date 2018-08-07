#' @title Test if an im object is a binary map
#' @export isbinarymap
#' @import spatstat
#' @param xiim an image object
#' @param requiretrue Logical. If TRUE then isbinarymap will error if xiim is NOT a binary map.
#' @description
#' Tests \code{xiim} to see if it is a binary map. The pixel values must be of logical type (TRUE, FALSE and NA only), or numerical with values of 0, 1 and NA.

#' @return Logical value. TRUE if \code{xiim} is a binary map. Otherwise FALSE.

#' @examples
#' #The following return TRUE
#' isbinarymap(as.im(heather$coarse, na.value = 0))
#' isbinarymap(as.im(heather$coarse, na.value = FALSE, value = TRUE))

#' #the following return FALSE
#' isbinarymap(as.im(heather$coarse, na.value = 0.2, value = 1))
#' isbinarymap(as.im(heather$coarse, na.value = 0, value = 1.5))

isbinarymap <- function(xiim, requiretrue = FALSE) {
    uvals <- unique(as.list(as.matrix(xiim)))
    if ( !all(  (uvals %in% c(0, 1)) | is.na(uvals))  && 
             !all((uvals %in% c(FALSE, TRUE, NA)) | is.na(uvals)) ) {
        isbinary <- FALSE
    } else {
        isbinary <- TRUE
    }
    if (requiretrue && !isbinary) {stop("Input xi has values other than 0, 1 or NA")}
    return(isbinary)
}
