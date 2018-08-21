context("Summarising Operations")

test_that("summary.imlist matches point wise summaries", {
  xy <- as.mask(square(1), eps = 0.1)
  unitname(xy) <- c("metre", "metres")
  zheights <- rnorm(100)
  zfun <- function(x, y){
    return(x * y)
  }
  ims <- lapply(zheights, function(z) {
           z * as.im(zfun, W = xy)
         })

  ims.summ <- summary.imlist(ims)
  expect_equivalent(ims.summ$mean[1, 1, drop = TRUE], mean(zheights * zfun(xy$xcol[[1]], xy$yrow[[1]])))
  expect_equivalent(ims.summ$var[5, 3, drop = TRUE], var(zheights * zfun(xy$xcol[[5]], xy$yrow[[3]])))
  expect_equivalent(ims.summ$max[5, 8, drop = TRUE], max(zheights * zfun(xy$xcol[[5]], xy$yrow[[8]])))
  expect_equivalent(ims.summ$min[3, 8, drop = TRUE], min(zheights * zfun(xy$xcol[[3]], xy$yrow[[8]])))

  expect_equal(unitname(ims[[1]]), unitname(xy))
  
  expect_error(ims.summ <- summary.imlist(ims, na.rm = TRUE))
})
