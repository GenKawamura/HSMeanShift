blurringMeanShiftOperator <- function( X, h=1, kernel="epanechnikovKernel" ){
	
	n.curves <- ncol( X )
	
	## compute distances
	distances <- as.matrix( dist( t( X ), diag=TRUE, upper=TRUE ) )
	
	## scale by bandwidth
	scaled.distances <- distances / h
	
	## evaluate kernel
	kernel <- get( kernel )
	kernel.values <- matrix( kernel( scaled.distances ), nrow=n.curves,
	ncol=n.curves ) 
	
	## weights denominators
	total.sum <- colSums( kernel.values )
	
	## weights
	kernel.weights <- kernel.values / total.sum

	## update
	new.X <- X%*%t( kernel.weights )
	
	output <- new.X
	
	return( new.X )
	
}


meanShiftOperator <- function( x, points, h=1,
kernel="epanechnikovKernel" ){
	
	## mean-shift operator
	
	## compute distances
	distances <- apply( points, 2, distanceFunction, y=x )
	
	## scale by bandwidth
	scaled.distances <- distances / h
	
	## evaluate kernel
	kernel <- get( kernel )
	kernel.values <- kernel( scaled.distances )
	
	## weights denominator
	total.sum <- sum( kernel.values )
	
	## mean-shift weights
	if( total.sum > 0 ){
		
		## update
		kernel.weights <- kernel.values / sum( kernel.values )
		output <- points %*% kernel.weights
		
	} else{
		
		output <- x
		
	}
	
	return( output )
	
}

###


gaussianKernel <- function( x ){
	
	## function to evaluate the asymmetric gaussian kernel	
	computeGaussianKernel <- function( y ){
	
		if( 0 <= y ){
		
			value <- 2 / 0.388 * dnorm( y / 0.388 )
		
		} else{
		
			value <- 0
		
		}
	
		return( value )
	
	}
	
	output <- sapply( x, computeGaussianKernel )
	
	return( output )
		
}


###

exponentialKernel <- function( x ){
	
	## function to evaluate the asymmetric exponential kernel	
	computeExponentialKernel <- function( y ){
	
		if( 0 <= y ){
		
			value <- dexp( y, rate=4.61 )
		
		} else{
		
			value <- 0
		
		}
	
		return( value )
	
	}
	
	output <- sapply( x, computeExponentialKernel )
	
	return( output )
		
}

###

cubicKernel <- function( x ){
	
	## function to evaluate the asymmetric cubic kernel	
	computeCubicKernel <- function( y ){
	
		if( 0 <= y && y<= 1 ){
		
			value <- 4 * ( 1 - y )^3
		
		} else{
		
			value <- 0
		
		}
	
		return( value )
	
	}
	
	output <- sapply( x, computeCubicKernel )
	
	return( output )
		
}

###

epanechnikovKernel <- function( x ){
	
	## function to evaluate the asymmetric Epanechnikov kernel	
	computeEpanechnikovKernel <- function( y ){
	
		if( 0 <= y && y<= 1 ){
		
			value <- 3 / 2 * ( 1 - y^2 )
		
		} else{
		
			value <- 0
		
		}
	
		return( value )
	
	}
	
	output <- sapply( x, computeEpanechnikovKernel )
	
	return( output )
		
}

###

distanceFunction <- function( x, y ){
	
	## function to compute the standard euclidean distance
	output <- sqrt( sum( ( x - y )^2 ) )
	
	return( output )
	
}

