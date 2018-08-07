context("Summary Functions")

test_that("summary.fvlist gives correct numbers for constant function", {
  xvals <- seq(0, 1, by = 0.1)
  yheights <- rnorm(100)
  fvlist <- lapply(yheights, function(y) as.fv(data.frame(xvals, y)) )
  fvlist.summ <- summary.fvlist(fvlist)
  expect_equal(fvlist.summ$maxy[[1]], max(yheights))
  expect_equal(fvlist.summ$miny[[1]], min(yheights))
  expect_equal(fvlist.summ$meany[[1]], mean(yheights))
  expect_equal(fvlist.summ$vary[[1]], var(yheights))
  expect_equal(fvlist.summ[1, fvnames(fvlist.summ, ".a"), drop = TRUE], 
    fvlist.summ[6, fvnames(fvlist.summ, ".a"), drop = TRUE]) 
})

#test when functions have multiple y columns
test_that("summary.fvlist works when multiple y columns around", {
  xvals <- seq(0, 1, by = 0.1)
  yheights <- rnorm(100)
  fvlist <- lapply(yheights, function(y) {
		   as.fv(data.frame(x = xvals, y = y, y2 = y^2, y3 = y^3))
    } )
  fvlist.summ <- summary.fvlist(fvlist)
  expect_length(fvnames(fvlist.summ, ".a"), 4 * length(fvnames(fvlist[[1]], ".a")))
  expect_equal(fvlist.summ$maxy[[1]], max(yheights))
  expect_equal(fvlist.summ$miny[[1]], min(yheights))
  expect_equal(fvlist.summ$meany[[1]], mean(yheights))
  expect_equal(fvlist.summ$vary[[1]], var(yheights))
  
  expect_equal(fvlist.summ$maxy2[[1]], max(yheights^2))
  expect_equal(fvlist.summ$miny2[[1]], min(yheights^2))
  expect_equal(fvlist.summ$meany2[[1]], mean(yheights^2))
  expect_equal(fvlist.summ$vary2[[1]], var(yheights^2))
  
  expect_equal(fvlist.summ$maxy3[[1]], max(yheights^3))
  expect_equal(fvlist.summ$miny3[[1]], min(yheights^3))
  expect_equal(fvlist.summ$meany3[[1]], mean(yheights^3))
  expect_equal(fvlist.summ$vary3[[1]], var(yheights^3))
})


#test when functions aren't constant
test_that("summary.fvlist operates when functions not constant", {
  xvals <- seq(0, 1, by = 0.1)
  yheights <- rnorm(100)
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
