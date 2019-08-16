#' @title Inner Product of Two Functions Represented as Images
#' @export innerprod.im
#'
#' @description Given two functions, f and g, that map from 2D space to 1D, and values of f and g represented as \code{im} objects
#'  then \code{innerprod.im} computes the (function space) inner product
#' \deqn{\int f(v) g(v) dv.}


#' @param A An \code{im} object containing function values representing function \eqn{f}.
#' @param B An \code{im} object containing function values representing function \eqn{g}.
#' @param outsideA The value of \eqn{f} outside the domain of \code{A}. Typically will be 0 or NA. Default is NA.
#' @param outsideB The value of \eqn{g} outside the domain of \code{B}. Typically will be \code{0} or \code{NA}. Default is \code{NA}.
#' @param na.rm Logical. If TRUE NA values are skipped in the summation of the product of the images.

#' @details This function harmonises the two input images, multiplies them together and returns the
#' integral of the resulting image.
#'  \code{outsideA} and \code{outsideB} are used to determine result if the inner product requires
#' values outside the domain of A or B. For example if \code{outsideA=0} and the
#' domain of \code{B} is larger than \code{A}'s domain then the inner product
#' can still be computed. However if \code{A} is \code{NA} outside (e.g. not
#' known/not provided) and the domain of \code{B} is larger than \code{A}'s
#' domain then the inner product can not be computed and the returned value is \code{NA}

#' @return 
#' If the inner product can be computed then returns sum(\code{A} * \code{B}), otherwise returns \code{NA}.

#' @examples
#' xi <- heather$coarse
#' covar <- plugincvc(xi, Frame(xi))
#' B <- setcov(square(1))
#' innerprod.im(covar, B, outsideB = 0)
innerprod.im <- function(A, B, outsideA = NA, outsideB = NA, na.rm = FALSE, method = "cubature"){
  ##assume that NA is unknown but finite so that mutliplying 0*NA gives 0
  integrationregion <- union.owin(Frame(A), Frame(B))
  #check if results will be NA due to outsideA and outsideB
  if (is.na(outsideA) && is.na(outsideB)){return(NA)}
  # if one is NA outside and the other is non-zero outside then result is NA
  if (is.na(outsideA) && (outsideB != 0)){return(NA)} 
  if (is.na(outsideB) && (outsideA != 0)){return(NA)} 
  #if one is non-zero and the other is non-zero then the result is Inf
  if ((!is.na(outsideA)) && (outsideA != 0) && (outsideB != 0)) {return(Inf)}
  if ((!is.na(outsideB)) && (outsideB != 0) && (outsideA != 0)) {return(Inf)}
 
  #outsideA is 0 or outsideB is 0 then ignore what happens outside the integration region 
  if ((!is.na(outsideA)) && (outsideA == 0)){integrationregion <- intersect.owin(integrationregion, Frame(A))}
  if ((!is.na(outsideB)) && (outsideB == 0)){integrationregion <- intersect.owin(integrationregion, Frame(B))}
  
  #now we have an integration region based on non-zero locations.
  #we have that at least one must be non-NA outside
  #we have that if one is NA and the other is non-zero outside then the result is NA
  #we that if one is non-zero and the other is non-zero then the result is Inf

  if ((method == "simple") || (requireNamespace("cubature") != TRUE)){
     intresult <- integration_trad(A, B, outsideA, outsideB, integrationregion)
  } else {
     intresult <- integration_cubature(A, B, outsideA, outsideB, integrationregion)$integral
  }
  return(intresult)
}


integration_trad <- function(A, B, outsideA, outsideB, integrationregion){
  harmgrid <- as.mask(integrationregion,
             eps = c(min(A$xstep, B$xstep), min(A$ystep, B$ystep)))
  A2 <- as.im(A, xy = harmgrid)
  A2[setminus.owin(integrationregion, Frame(A))] <- outsideA
  B2 <- as.im(B, xy = harmgrid)
  B2[setminus.owin(integrationregion, Frame(B))] <- outsideB
  prdimg <- eval.im(A2 * B2, harmonize = FALSE)
  return(sum(prdimg[, ], na.rm = na.rm) * prdimg$xstep * prdimg$ystep)
}

integration_cubature <- function(A, B, outsideA, outsideB, integrationregion, tol = 1E-3){
  if (requireNamespace("cubature") != TRUE){
     stop("Cubature package must be installed to integrate using integration_cubature")
  }
  
  tmpfunA <- as.function.im(A)
  tmpfunB <- as.function.im(B)
  #vectorised cubature functions needs function that take an matrix with each point a column
  integrand <- function(arg){
    insideA <- inside.owin(x = arg[1, ], y = arg[2, ], w = Window(A))
    insideB <- inside.owin(x = arg[1, ], y = arg[2, ], w = Window(B))
    outA_tmp <- matrix(tmpfunA(arg[1, insideA], arg[2, insideA]), ncol = ncol(arg)) 
    outA <- matrix(outsideA, nrow = nrow(outA_tmp), ncol = ncol(arg), byrow = FALSE)
    outA[, insideA] <- outA_tmp

    outB_tmp <- matrix(tmpfunB(arg[1, insideB], arg[2, insideB]), ncol = ncol(arg)) 
    outB <- matrix(outsideB, nrow = nrow(outB_tmp), ncol = ncol(arg), byrow = FALSE)
    outB[, insideB] <- outB_tmp

    return(outA * outB) 
  }
  out <- cubature::hcubature(f = integrand, 
                   lowerLimit = c(integrationregion$xrange[[1]], integrationregion$yrange[[1]]),
                   upperLimit = c(integrationregion$yrange[[2]], integrationregion$yrange[[2]]),
                   tol = tol,  #stops when integral accurate to tol * integral value
                   vectorInterface = TRUE)
  return(out)
}
