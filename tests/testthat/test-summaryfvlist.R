context("Summary Functions")

test_that("summary.fvlist gives correct numbers for constant function", {
  xvals <- seq(0, 1, by = 0.1)
  yheights <- rnorm(100, mean = 5)
  fvlist <- lapply(yheights, function(y) as.fv(data.frame(xvals, y)) )
  fvlist.summ <- summary.fvlist(fvlist)
  expect_equal(fvlist.summ$maxy[[1]], max(yheights))
  expect_equal(fvlist.summ$miny[[1]], min(yheights))
  expect_equal(fvlist.summ$meany[[1]], mean(yheights))
  expect_equal(fvlist.summ$vary[[1]], var(yheights))
  expect_equal(fvlist.summ[1, fvnames(fvlist.summ, ".a"), drop = TRUE], 
    fvlist.summ[6, fvnames(fvlist.summ, ".a"), drop = TRUE]) 
  expect_equal(fvlist.summ$nay[[1]], 0)
})

#test when functions have multiple y columns
test_that("summary.fvlist works when multiple y columns around", {
  xvals <- seq(0, 1, by = 0.1)
  yheights <- rnorm(100)
  fvlist <- lapply(yheights, function(y) {
		   as.fv(data.frame(x = xvals, y = y, y2 = y^2, y3 = y^3))
    } )
  fvlist.summ <- summary.fvlist(fvlist)
  expect_length(fvnames(fvlist.summ, ".a"), 5 * length(fvnames(fvlist[[1]], ".a")))
  expect_equal(fvlist.summ$maxy[[1]], max(yheights))
  expect_equal(fvlist.summ$miny[[1]], min(yheights))
  expect_equal(fvlist.summ$meany[[1]], mean(yheights))
  expect_equal(fvlist.summ$vary[[1]], var(yheights))
  expect_equal(fvlist.summ$nay[[1]], 0)
  
  expect_equal(fvlist.summ$maxy2[[1]], max(yheights^2))
  expect_equal(fvlist.summ$miny2[[1]], min(yheights^2))
  expect_equal(fvlist.summ$meany2[[1]], mean(yheights^2))
  expect_equal(fvlist.summ$vary2[[1]], var(yheights^2))
  expect_equal(fvlist.summ$nay2[[1]], 0 )
  
  expect_equal(fvlist.summ$maxy3[[1]], max(yheights^3))
  expect_equal(fvlist.summ$miny3[[1]], min(yheights^3))
  expect_equal(fvlist.summ$meany3[[1]], mean(yheights^3))
  expect_equal(fvlist.summ$vary3[[1]], var(yheights^3))
  expect_equal(fvlist.summ$nay3[[1]], 0)
})


#test when functions aren't constant
test_that("summary.fvlist operates when functions not constant", {
  xvals <- seq(0, 1, by = 0.1)
  yheights <- rnorm(100, mean = 5)
  fvlist <- lapply(yheights, function(y) {
		   as.fv(data.frame(x = xvals, y = y * xvals, y2 = (y * xvals)^2))} )
  fvlist.summ <- summary.fvlist(fvlist) 
  
  expect_equal(fvlist.summ$maxy[[5]], max(yheights * xvals[[5]]))
  expect_equal(fvlist.summ$miny[[5]], min(yheights * xvals[[5]]))
  expect_equal(fvlist.summ$meany[[5]], mean(yheights * xvals[[5]]))
  expect_equal(fvlist.summ$vary[[5]], var(yheights * xvals[[5]]))

  expect_equal(fvlist.summ$maxy2[[5]], max((yheights * xvals[[5]])^2))
  expect_equal(fvlist.summ$miny2[[5]], min((yheights * xvals[[5]])^2))
  expect_equal(fvlist.summ$meany2[[5]], mean((yheights * xvals[[5]])^2))
  expect_equal(fvlist.summ$vary2[[5]], var((yheights * xvals[[5]])^2))
})

#test when some functions are NA for some parts of the abscissa
test_that("summary.fvlist operates when functions contain NA values", {
  xvals <- seq(0, 1, by = 0.1)
  yheights_a <- rnorm(3, mean = 2)
  yheights_b <- rnorm(3, mean = 5)
  yheights <- c(yheights_a, yheights_b)
  fvlist_a <- lapply(yheights_a, function(y) {
    as.fv(data.frame(x = xvals, y = y * xvals, y2 = (y * xvals)^2))} )
  fvlist_b <- lapply(yheights_b, function(y) {
    df <- data.frame(x = xvals, y = y * xvals, y2 = (y * xvals)^2)
    df$y[xvals > 0.5] <- NA
    df$y2[xvals < 0.25 | xvals > 0.75] <- NA
    as.fv(df)
  } )
  fvlist <- c(fvlist_a, fvlist_b)
  
  #na.rm = FALSE
  fvlist.summ <- summary.fvlist(fvlist, na.rm = FALSE)
  expect_true(is.na(fvlist.summ$meany[[7]]))
  expect_false(is.na(fvlist.summ$meany2[[7]]))
  
  #na.rm = TRUE
  fvlist.summ <- summary.fvlist(fvlist, na.rm = TRUE)
  
  expect_equal(fvlist.summ$maxy[[2]], max(yheights * xvals[[2]]))
  expect_equal(fvlist.summ$miny[[2]], min(yheights * xvals[[2]]))
  expect_equal(fvlist.summ$meany[[2]], mean(yheights * xvals[[2]]))
  expect_equal(fvlist.summ$vary[[2]], var(yheights * xvals[[2]]))
  expect_equal(fvlist.summ$nay[[2]], 0)
  
  expect_equal(fvlist.summ$maxy[[8]], max(yheights_a * xvals[[8]]))
  expect_equal(fvlist.summ$miny[[8]], min(yheights_a * xvals[[8]]))
  expect_equal(fvlist.summ$meany[[8]], mean(yheights_a * xvals[[8]]))
  expect_equal(fvlist.summ$vary[[8]], var(yheights_a * xvals[[8]]))
  expect_equal(fvlist.summ$nay[[8]], length(fvlist_b))
  
  expect_equal(fvlist.summ$maxy2[[9]], max((yheights_a * xvals[[9]])^2))
  expect_equal(fvlist.summ$miny2[[9]], min((yheights_a * xvals[[9]])^2))
  expect_equal(fvlist.summ$meany2[[9]], mean((yheights_a * xvals[[9]])^2))
  expect_equal(fvlist.summ$vary2[[9]], var((yheights_a * xvals[[9]])^2))
  expect_equal(fvlist.summ$nay2[[9]], length(fvlist_b))
})


