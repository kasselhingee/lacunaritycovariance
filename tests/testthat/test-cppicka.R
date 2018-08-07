context("cppicka")

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

