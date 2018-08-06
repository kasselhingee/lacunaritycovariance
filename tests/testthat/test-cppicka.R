context("cppicka")

#cppikca is already tested indirectly through covar().
# it would be very cool in the future to use my analytical computations for rectangles to test cppicka
cpest <- cppicka(heather$coarse, obswin = Frame(heather$coarse))
phat <- coverageprob(heather$coarse, obswin = Frame(heather$coarse))

test_that("cppicka output is unsymmetric", {
  diffsymm <- cpest - reflect.im(cpest)
  diffsymm[setcov(Frame(heather$coarse)) > 0.2 *area.owin(Frame(heather$coarse))] <- NA #ignore small vectors
  expect_gt(mean(abs(diffsymm)), 0.01)
})

test_that("cppicka matchs phat at centre", {
  cppickacentre <- cpest[ ppp(x = 0, y = 0, W = disc(r = 1))]
  expect_equal(cppickacentre, phat)
})
