context("test-rcpp_functions")

test_that("Testing Rcpp functions if it has the same functionality as R functions have", {

  ## Check distanceFunction
  x <- runif(100, 0, 10)
  y <- runif(100, 0, 10)
  
  expect_equal(distanceFunction(x,x), 0)
  expect_equal(RcppDistanceFunction(x,x), 0)
  expect_equal(distanceFunction(x,y), RcppDistanceFunction(x,y))

  ## Check EpanechnikovKernel
  x.1 <- runif(1000, 0, 2)
  y.1 <- runif(1000, 0, 2)
  expect_equal(epanechnikovKernel(x.1), RcppEpanechnikovKernel(x.1))
  expect_equal(epanechnikovKernel(y.1), RcppEpanechnikovKernel(y.1))

})

