context("test-rcpp_epanechnikov_kernel")

test_that("Testing Rcpp epanechnikovKernel function", {
  ## Setting up
  source("setup.R")
  
  ##----------------------------------
  ## Kernel
  ##----------------------------------
  kernel <- "epanechnikovKernel"

  
  ## Check
  res1 <- meanShiftOperator(x, points, h, kernel)
  res2 <- RcppMeanShiftOperator(x, points, h, kernel)
  expect_equal(sum(res1), sum(res2))
  

  ## Check
  res3 <- RcppMeanShiftAlgorithm(x, points, h, kernel, tol.stop)
  res4 <- meanShiftAlgorithm(x, points, h, kernel, tol.stop)

  expect_equal(sum(res3), sum(res4))
  
  ## Running msClustering
  rcpp.clustering <- msClustering( iris.data, h, kernel=kernel, multi.core=TRUE, legacy.mode=FALSE )
  non.rcpp.clustering <- msClustering( iris.data, h, kernel=kernel, multi.core=TRUE, legacy.mode=TRUE )
  
  expect_equal(sum(rcpp.clustering$components), sum(non.rcpp.clustering$components))
})

