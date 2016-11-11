## an example using the iris dataset
## help( iris )

## prepare data matrix (a subset of the iris dataset)
set.seed( 2 )
indices <- sample( 1:nrow( iris ), 80 )
iris.data <- t( iris[indices,c( "Sepal.Length", "Sepal.Width" )] )

## run mean shift algorithm
clustering <- msClustering( iris.data, h=0.8 )
print( clustering )

## plot the clusters
## Not run:
plot( iris.data[1,], iris.data[2,], col=clustering$labels+2, cex=0.8,
     pch=16, xlab="Sepal.Length", ylab="Sepal.Width" )
points( clustering$components[1,], clustering$components[2,],
       col=2+( 1:ncol( clustering$components ) ), cex=1.8, pch=16 )
## End(Not run)

## using multiple cores (2)
## Not run:
options( mc.cores=2 )
clustering.mc <- msClustering( iris.data, multi.core=TRUE )
## End(Not run)
