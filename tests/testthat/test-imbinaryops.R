context("Summarising Operations")

test_that("im binary operations are performing correctly", {
  win <- as.im(square(10), eps = 0.25)
  #create two ims with unequal pixel arrays
  ims <- solist(as.im(square(1), W = win, xy = win, na.replace = 0),
                as.im(shift.owin(square(1), vec = c(0.75, 0.75)), W = win[square(4), drop = TRUE],
                      xy = win[square(4), drop = TRUE], na.replace = 0))
  expect_error(Pmax.im(ims[[1]], ims[[2]]), regexp = "not compatible")
  expect_error(Pmin.im(ims[[1]], ims[[2]]), regexp = "not compatible")
  
  ims <- as.solist(do.call(harmonize.im, args = ims))

  expect_equal(sum(ims[[1]]) + sum(ims[[2]]) - sum(Pmax.im(ims[[1]], ims[[2]])), sum(Pmin.im(ims[[1]], ims[[2]])))
  
  expect_equal(sum(Add.im(ims[[1]], ims[[2]])), sum(ims[[1]]) + sum(ims[[2]]))
  expect_equal(sum(Square.im(Add.im(ims[[1]], ims[[2]]))), 34)
})

