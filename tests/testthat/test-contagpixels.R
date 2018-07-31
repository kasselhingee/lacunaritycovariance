context("Contagion Computations")

test_that("Pixel Adjacency Contagion gives correct results on toy images", {
#should have 0 nbrs of Xi and Xi, and 4 nbrs of xi and not xi
#should have 4 nbrs of not xi and xi, and lots of not xi and not xi
testmat <- matrix(FALSE,nrow=5,ncol=5)
testmat[2,2] <- TRUE
xi <- owin(mask = testmat)
w <- Frame(xi)
adjmat <- adjacency(xi, Frame(xi))
expect_equal(adjmat,
	     matrix(c(0, 4, 4, 72), nrow = 2, ncol = 2, byrow = TRUE),
	     check.attributes = FALSE)

#should have 2 nbrs of Xi and Xi, and 6 nbrs of xi and not xi
#should have 6 nbrs of not xi and xi, and lots of not xi and not xi
testmat <- matrix(FALSE,nrow=5,ncol=5)
testmat[2,2] <- TRUE
testmat[3,2] <- TRUE
xi <- owin(mask = testmat)
w <- Frame(xi)
adjmat <- adjacency(xi, Frame(xi))
expect_equal(adjmat,
	     matrix(c(2, 6, 6, 66), nrow = 2, ncol = 2, byrow = TRUE),
	     check.attributes = FALSE)

#should have 2 nbrs of Xic and Xic, and 6 nbrs of xic and xi
#should have 6 nbrs of xi and xic, and lots of xi and xi
testmat <- matrix(TRUE,nrow=5,ncol=5)
testmat[2,2] <- FALSE
testmat[2,3] <- FALSE
xi <- owin(mask = testmat)
w <- Frame(xi)
adjmat <- adjacency(xi, Frame(xi))
expect_equal(adjmat,
	     matrix(c(66, 6 , 6, 2), nrow = 2, ncol = 2, byrow = TRUE),
	     check.attributes = FALSE)


#minimum contagion by making every pixel have half nbrs in Xi and half outside xi
vectmp <- vector(mode="logical",length=100)
vectmp[(1:length(vectmp) %%2)==0] <- TRUE
vectmp[(1:length(vectmp) %%2)!=0] <- FALSE
testmat <- matrix(vectmp,nrow=10,ncol=10)
xi <- owin(mask = testmat)
expect_equal(contagpixelgrid(xi,Frame(xi),normalise=TRUE), 0)

#maximum contagion
vectmp <- vector(mode="logical",length=10000)
vectmp[3:length(vectmp)] <- TRUE
vectmp[1:2] <- FALSE #two points to avoid zero adjacency
testmat <- matrix(vectmp,nrow=100,ncol=100)
xi <- owin(mask = testmat)
expect_equal(contagpixelgrid(xi,Frame(xi),normalise=TRUE),
	     1, tolerance = 0.01)
#result should be close to 1
})