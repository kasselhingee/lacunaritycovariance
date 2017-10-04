context("RACS simulation")

test_that("Boolean model with log normal disc has area fraction consistent with past", {
  ##use seed to fix the simulation
  w <- owin(xrange = c(0, 10), yrange = c(0, 10))
  xi <- rblnd(w, 2, 1, -1, 0.5, seed = 36)
  expect_equal(area.owin(xi), 49.2101, tolerance = 1E-7)
})
