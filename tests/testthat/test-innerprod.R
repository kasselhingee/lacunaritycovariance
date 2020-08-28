context("Internals")

test_that("Inner product gives correct error results",  {
  A <- as.im(function(x,y) {sin(x/2.0)},W=square(4*pi),eps=0.5)
  B <- as.im(function(x,y) {sin(x)},W=square(2*pi),eps=0.5)
  expect_true(is.na(innerprod.im(A, B, outsideA = NA, outsideB = NA)))
  expect_false(is.na(innerprod.im(A, B, outsideA = 0, outsideB = 0)))
  
  expect_true(is.na(innerprod.im(A, B, outsideA = 0, outsideB = NA)))
  
  expect_false(is.na(innerprod.im(A, B, outsideA = NA, outsideB = 0)))
  
  expect_true(is.na(innerprod.im(A, B, outsideA = 1, outsideB = NA)))
  expect_true(is.na(innerprod.im(A, B, outsideA = NA, outsideB = 1)))
  expect_false(is.na(innerprod.im(A, B, outsideA = 0, outsideB = 1)))
  expect_true(is.infinite(innerprod.im(A, B, outsideA = 1, outsideB = 1)))
  
  A <- as.im(function(x,y) {sin(x/2.0)},W=square(4*pi),eps=0.01)
  B <- as.im(function(x,y) {sin(x)},W=square(2*pi),eps=0.01)
  #innerproduct should be close to 0 (orthogonal):
  #innerprod.im(as.im(function(x,y) {sin(x)},W=square(7*pi),eps=0.01),as.im(function(x,y) {sin(2*x)},W=square(2*pi),eps=0.01))
  expect_lt(abs(innerprod.im(A, B, outsideB = 0)), 0.1)
  expect_gt(abs(innerprod.im(A, 1 + 0 * B, outsideB = 0)), 0.1)
  expect_lt(abs(innerprod.im(A, 1 + 0 * B, outsideA = 0, outsideB = 1)), 0.1)
 
  skip_on_cran() #mostly a repeat of above tests 
  A <- as.im(function(x,y) {sin(x/2.0)}, W = owin(xrange = c(0, 4*pi), yrange = c(0, pi)), eps = 0.5)
  B <- as.im(function(x,y) {sin(x)}, W = square(2 * pi), eps = 0.5)
  expect_true(is.na(innerprod.im(A, B, outsideA = NA, outsideB = 0, na.replace = FALSE)))
  expect_true(is.na(innerprod.im(A, B, outsideB = NA, outsideA = 0, na.replace = FALSE)))
  expect_true(is.infinite(innerprod.im(A, B, outsideB = 0.1, outsideA = 0.1)))
  expect_true(is.na(innerprod.im(A, B, outsideB = 0.1, outsideA = NA)))
  
  A <- as.im(function(x,y) {sin(x/2.0)}, W = owin(xrange = c(0, 4*pi), yrange = c(0, pi)), eps = 0.01)
  B <- as.im(function(x,y) {sin(x)}, W = square(2 * pi), eps = 0.01)
  expect_lt(abs(innerprod.im(A, B, outsideA = 0, outsideB = 0)), 0.1)
  expect_gt(abs(innerprod.im(A, 1 + 0 * B, outsideA = 0, outsideB = 0)), 0.1)
  expect_lt(abs(innerprod.im(1 + 0 * A, B, outsideA = 0, outsideB = 0)), 0.1)
})
