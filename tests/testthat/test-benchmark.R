context("test-benchmark")

test_that("Benchmarking Rcpp function", {

  ## Setting up
  options( mc.cores=8 )
  
  set.seed( 2 )
  indices <- sample( 1:nrow( iris ), 150 )
  iris.data <- t( iris[indices,c( "Sepal.Length", "Sepal.Width" )] )

  points <- as.matrix(iris.data)
  x <- lapply( apply( points, 2, list ), unlist )[[1]]

  kernel <- "epanechnikovKernel"
  h <- 0.8
  tol.stop <- 1e-6

  
  ## Benchmarking msClustering
  require(rbenchmark)
  a <- benchmark(msClustering( iris.data, h, multi.core=FALSE, legacy.mode=TRUE ),
            msClustering( iris.data, h, multi.core=FALSE, legacy.mode=FALSE ),
            columns=c("test", "replications", "elapsed", "relative"),
            order="relative",
            replications=10)
  print(a)

  b <- benchmark(msClustering( iris.data, h, multi.core=TRUE, legacy.mode=TRUE ),
            msClustering( iris.data, h, multi.core=TRUE, legacy.mode=FALSE ),
            columns=c("test", "replications", "elapsed", "relative"),
            order="relative",
            replications=10)
  print(b)
})

