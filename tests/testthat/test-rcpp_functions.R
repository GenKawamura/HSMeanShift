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

  ## Setting up
  options( mc.cores=8 )
  
  set.seed( 2 )
  indices <- sample( 1:nrow( iris ), 150 )
  iris.data <- t( iris[indices,c( "Sepal.Length", "Sepal.Width" )] )

  points <- as.matrix(iris.data)
  x <- lapply( apply( points, 2, list ), unlist )[[1]]


  ##----------------------------------
  ## Cubic Kernel
  ##----------------------------------

  kernel <- "cubicKernel"
  h <- 0.8
  tol.stop <- 1e-6

  
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


  ##----------------------------------
  ## Epanechnikov Kernel
  ##----------------------------------
  
  kernel <- "epanechnikovKernel"
  h <- 0.8
  tol.stop <- 1e-6

  
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

