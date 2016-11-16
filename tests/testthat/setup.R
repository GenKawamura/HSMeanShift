options( mc.cores=8 )

set.seed( 2 )
indices <- sample( 1:nrow( iris ), 150 )
iris.data <- t( iris[indices,c( "Sepal.Length", "Sepal.Width" )] )

points <- as.matrix(iris.data)
x <- lapply( apply( points, 2, list ), unlist )[[1]]

h <- 0.8
tol.stop <- 1e-6
