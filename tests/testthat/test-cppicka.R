context("Estimation - Covariance")

#cppikca is already tested indirectly through covar().
# it would be very cool in the future to use my analytical computations for rectangles to test cppicka

#creating a test that involves NA values in an image
w <- setminus.owin(Frame(heather$coarse), square(5)) 
xiowin <- intersect.owin(heather$coarse, w)
cpest <- cppicka(xiowin, obswin = w)
phat <- coverageprob(xiowin, obswin = w)


test_that("cppicka output is unsymmetric", {
  diffsymm <- cpest - reflect.im(cpest)
  diffsymm[setcov( w ) > 0.2 *area.owin( w)] <- NA #ignore small vectors
  expect_gt(mean(abs(diffsymm)), 0.01)
})

test_that("cppicka matchs phat at centre", {
  cppickacentre <- cpest[ ppp(x = 0, y = 0, W = disc(r = 1))]
  expect_equal(cppickacentre, phat)
})

test_that("cppicka gives expected results for a toy image", {
  #test that it is doing what we expect for a single vector
  xi <- shift.owin(square(r = 1), vec = c(2, 2))
  win <- square(r = 4)
  cpp1 <- cppicka(xi, obswin = win)
  #expect cppicka at v = c(2,2) to return 1/(2*2)
  expect_equal(cpp1[ppp(x = 2, y = 2, window = Frame(cpp1))],
               1/4  )
  #expect cppicka at v = c(-2, -2) to return 0
  expect_equal(cpp1[ppp(x = -2, y = -2, window = Frame(cpp1))],
               0  )
  #expect cppicka at v = c(2, 0) to return 1 / 8
  expect_equal(cpp1[ppp(x = 2, y = 0, window = Frame(cpp1))],
               1 / 8  )
})
