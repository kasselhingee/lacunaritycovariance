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
#' @param na.replace Logical. If TRUE NA values in \code{A} and \code{B} are replaced by \code{outsideA} and \code{outsideB}, respectively. This allows the cubature integration to be performed (roughly) twice as quickly.
#' @param method Either "cubature" or "harmonisesum". The former uses [cubature::cubintegrate()], 
#' the latter harmonises the images using [spatstat::as.mask()] and sums the product of \code{A} and \code{B}.

#' @details This function uses the package \pkg{cubature} to integrate the multiple of the two images. If \pkg{cubature} is not available then it harmonises the two input images, multiplies them together and returns the
#' integral of the resulting image.
#'  \code{outsideA} and \code{outsideB} are used to determine result if the inner product requires
#' values outside the domain of A or B. For example if \code{outsideA=0} and the
#' domain of \code{B} is larger than \code{A}'s domain then the inner product
#' can still be computed. However if \code{A} is \code{NA} outside (e.g. not
#' known/not provided) and the domain of \code{B} is larger than \code{A}'s
#' domain then the inner product can not be computed and the returned value is \code{NA}
#'
#' The \code{harmonisesum} method appears to perform well when used by [gblcc()], but produces numerical artefacts when used by [gblc()] and [gblg()]. The \code{cubature} method takes longer to compute and is more accurate for functions (A or B) that are non-zero for large vectors. This makes it suitable for use by [gblc()] and [gblg()].

#' @return 
#' If the inner product can be computed then returns sum(\code{A} * \code{B}), otherwise returns \code{NA}.

#' @examples
#' xi <- heather$coarse
#' covar <- plugincvc(xi, Frame(xi))
#' B <- setcov(square(1))
#' innerprod.im(covar, B, outsideB = 0)
innerprod.im <- function(A, B, outsideA = NA, outsideB = NA, na.replace = TRUE, method = "cubature"){
  ##assume that NA is unknown but finite so that mutliplying 0*NA gives 0
  integrationrectangle <- union.owin(Frame(A), Frame(B))
  #check if results will be NA due to outsideA and outsideB
  if (is.na(outsideA) && is.na(outsideB)){return(NA)}
  # if one is NA outside and the other is non-zero outside then result is NA
  if (is.na(outsideA) && (outsideB != 0)){return(NA)} 
  if (is.na(outsideB) && (outsideA != 0)){return(NA)} 
  #if one is non-zero and the other is non-zero then the result is Inf
  if ((!is.na(outsideA)) && (outsideA != 0) && (outsideB != 0)) {return(Inf)}
  if ((!is.na(outsideB)) && (outsideB != 0) && (outsideA != 0)) {return(Inf)}
 
  #outsideA is 0 or outsideB is 0 then ignore what happens outside the integration region 
  if ((!is.na(outsideA)) && (outsideA == 0)){integrationrectangle <- intersect.owin(integrationrectangle, Frame(A))}
  if ((!is.na(outsideB)) && (outsideB == 0)){integrationrectangle <- intersect.owin(integrationrectangle, Frame(B))}
  
  #now we have an integration region based on non-zero locations.
  #we have that outside the region must be zero (assuming 0 * NA = 0)
  #Other cases have already been computed and the function exits:
  # (a) if one is NA and the other is non-zero outside, then the result is NA
  # (b) if both are non-zero (and non-NA) outside then the result is Inf

  # the images may still have NA values inside their support windows, and there is still the choice between the different integration methods.
  if (method == "cubature"){
    intresult <- innerprod.cub(A, B, outsideA, outsideB, integrationrectangle, na.replace = na.replace)$integral
  } else if (method == "harmonisesum") {
    intresult <- integration_trad(A, B, outsideA, outsideB, integrationrectangle, na.replace = na.replace)
  } else {
    stop(paste("Method", method, "did not match available methods."))
  }
  return(intresult)
}


integration_trad <- function(A, B, outsideA, outsideB, integrationrectangle, na.replace){
  if (na.replace){
    A[is.na(as.matrix(A))] <- outsideA
    B[is.na(as.matrix(B))] <- outsideB
  }
  harmgrid <- as.mask(integrationrectangle,
             eps = c(min(A$xstep, B$xstep), min(A$ystep, B$ystep)))
  A2 <- as.im(A, xy = harmgrid)
  A2[setminus.owin(integrationrectangle, Frame(A))] <- outsideA
  B2 <- as.im(B, xy = harmgrid)
  B2[setminus.owin(integrationrectangle, Frame(B))] <- outsideB
  prdimg <- eval.im(A2 * B2, harmonize = FALSE)
  return(sum(prdimg[, ]) * prdimg$xstep * prdimg$ystep)
}

# dealing correctly with NA values (na.replace = FALSE) can double the computation time.
innerprod.cub <- function(A, B, outsideA, outsideB, integrationrectangle, na.replace, tol = 1E-3){
  if (requireNamespace("cubature") != TRUE){
     stop("Cubature package must be installed to integrate using method cubature")
  }
  
  if (na.replace){
  #vectorised cubature functions needs function that take an matrix with each point a column
    integrand <- function(arg){
      outA <- matrix(A[list(x = arg[1, ], y = arg[2, ]), drop = FALSE], ncol = ncol(arg))
      outA <- replace(outA, is.na(outA), outsideA)
      outB <- matrix(B[list(x = arg[1, ], y = arg[2, ]), drop = FALSE], ncol = ncol(arg))
      outB <- replace(outB, is.na(outB), outsideB)
      return(outA * outB)
    }
  } else if (!na.replace){
    tmpfunA <- as.function.im(A)
    tmpfunB <- as.function.im(B)
    #vectorised cubature functions needs function that take an matrix with each point a column
    integrand <- function(arg){
      insideA <- inside.owin(x = arg[1, ], y = arg[2, ], w = Window(A))
      insideB <- inside.owin(x = arg[1, ], y = arg[2, ], w = Window(B))
      outA <- matrix(outsideA, nrow = 1, ncol = ncol(arg), byrow = FALSE)
      if (sum(insideA) > 0){
        outA[, insideA] <- matrix(tmpfunA(arg[1, insideA], arg[2, insideA]), ncol = sum(insideA)) 
      }
      outB <- matrix(outsideB, nrow = 1, ncol = ncol(arg), byrow = FALSE)
      if (sum(insideB) > 0){
        outB[, insideB] <- matrix(tmpfunB(arg[1, insideB], arg[2, insideB]), ncol = sum(insideB)) 
      }
      return(outA * outB) 
    }
  }
  out <- cubature::cubintegrate(f = integrand, 
                   lower = c(integrationrectangle$xrange[[1]], integrationrectangle$yrange[[1]]),
                   upper = c(integrationrectangle$xrange[[2]], integrationrectangle$yrange[[2]]),
                   relTol = tol,  #stops when integral accurate to tol * integral value
                   method = "pcubature",
                   nVec = 1024)
  return(out)
}
