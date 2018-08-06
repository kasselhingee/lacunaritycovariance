context("Median ISE")

test_that("Median ISE gives the median value and the expected value", {
  rvals <- seq(0, 10, by = 0.01)
  yvals <- rvals
  fvlist <- lapply(c(0.1, 0.8, 1, 1.2, 1.4),
		   function(x) 
		   as.fv(data.frame(r = rvals, 
			    y = x * yvals)))
  reffv <- as.fv(data.frame(r = rvals, y = 0 * rvals))
  medise <- median_ise.fvlist(fvlist, reffv, domainlim = c(0, 10),
		    fixeddomain = TRUE,
		    avoverdomain = FALSE)
  expect_equal(medise, 10^3 / 3, tol = 1E-3)
})
