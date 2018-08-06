context("Kernel Smooth")

test_that("kernel smooth is believable for a step function", {
  stepim <- im(matrix(0,nrow=1000,ncol=1000),xrange=c(-1,1),yrange=c(-1,1))
  stepim[owin(xrange=c(-0.5,0.5),yrange=c(-0.5,0.5))] <- 1
  smstepim <- kernelsmooth(stepim, 0.5)
  #since Epanechnikov kernel bandwidth means support of kernel finishes at distances of 0.5
  samplevals <- smstepim[ppp(x = c(-1, -0.5, 0), y = c(0, 0, 0), window = Frame(stepim))]
  expect_equal(samplevals, c(0, 0.5, 1), tol = 1E-2)

  #integral should stay same since stepim is 0 within 0.5 of the boundary
  expect_equal(integral(stepim), integral(smstepim), tol = 1E-6)
})
