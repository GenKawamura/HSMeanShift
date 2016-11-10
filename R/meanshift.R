meanShiftAlgorithm <- function( x, points, h=1, kernel="epanechnikovKernel",
tol.stop=1e-6 , legacy.mode=FALSE){
  if (legacy.mode){
    close.enough <- FALSE
    
    old.x <- x
    
    ## while the trajectory has not converged
    ## (update produced a shift larger than 'tol.stop')
    while( !close.enough ) {
      
      ## apply mean-shift operator and update
      new.x <- meanShiftOperator( x=old.x, points=points, h=h,
                                 kernel=kernel)
      
      distance <- distanceFunction( old.x, new.x )
      
      old.x <- new.x
      
      close.enough <- ( distance < tol.stop )
      
    }

    return( new.x )

  } else {
    ## Rcpp implementation of meanShiftAlgorithm.
    ## This method should be much faster than ordinary R implementation
    output <- RcppMeanShiftAlgorithm(x, points, h, kernel, tol.stop)
    rownames(output) <- rownames(points)
    return(output)
  }
}

###

meanShiftAlgorithmAll <- function( X, h=NULL, kernel="epanechnikovKernel",
tol.stop=1e-6, multi.core=FALSE , legacy.mode=FALSE){
	
	if( is.null( h ) ){
		
		h <- quantile( dist( t( X ) ), 0.3 )
		
	}
	
	if( multi.core ){
		
		## MULTICORE REQUIRES 'parallel' LIBRARY
				
		X.list <- lapply( apply( X, 2, list ), unlist )
		
		multi.core.output <- mclapply( X.list, meanShiftAlgorithm,
		points=X, h=h, kernel=kernel, tol.stop=tol.stop, legacy.mode=legacy.mode )
		
		output <- do.call( cbind, multi.core.output )
                
	} else{
		
		M <- X
		n <- ncol( X )
		
		pb <- txtProgressBar( min=0, max=n, style=3 )
		
		for( i in 1:n ){
			
			M[,i] <- meanShiftAlgorithm( x=X[,i], points=X, h=h,
			kernel=kernel, tol.stop=tol.stop, legacy.mode=legacy.mode)
			
			setTxtProgressBar( pb, i )
			
		}
		
		close( pb )
		
		output <- M
	
	}

	message( "\nMean-shift algorithm ran successfully.\n")
	
	return( output )
	
}


#' Function to perform clustering using the mean shift algorithm.
#'
#' This function implements the mean shift algorithm. The algorithm locates the modes of a kernel density estimator and associates each data point to exactly one of the modes, thus effectively clustering the data.
#'
#' @useDynLib MeanShift
#' @export
#' @param X
#' @param h
#' @param kernel
#' @param tol.stop
#' @param tol.epsilon
#' @param multi.core
#' @param rcpp.mode
#' @return
msClustering <- function( X, h=NULL, kernel="epanechnikovKernel",
tol.stop=1e-6, tol.epsilon=1e-3, multi.core=FALSE , legacy.mode=FALSE){
	
	# minimal input checking
	X <- as.matrix( X )
		
	if( ncol( X ) <= 1 ){
		
		message( "The input matrix X has only one column: ",
		"returning input.")
		return( X )
	}
	
	if( !is.element( kernel, paste( c( "epanechnikov", "cubic", 
	"gaussian", "exponential"), "Kernel", sep="" ) ) ){
		
		stop( "Invalid kernel name.")
		
	}
	
	if( !is.null( h ) && h <= 0 ){
		
		stop( "The bandwidth must be strictly positive." )
		
	}
		
	if( tol.stop <= 0 || tol.epsilon <= 0 ){
		
		stop( "All tolerances must be strictly positive.")
		
	}
		
	## run mean-shift algorithm
	message( "\nRunning mean-shift algorithm...\n" )
	
	if( multi.core ){
		
		n.cores <- getOption( "mc.cores" )
		
		if( is.null( n.cores ) ){
			
			readInteger <- function(){
				
				n <- readline( "Enter number of cores: " )
				
				n <- as.integer( n )
				
				return( n )
				
			}
			
			n.cores <- readInteger()
			
			if( n.cores < 1 ){
				
				n.cores <- 1
				options( mc.cores=n.cores )
				
				cat( "\n" )
								
				warning( "\nInvalid choice for the number ",
				"of cores: 'mc.cores' option set to 1. To change, ",
				"use options( mc.cores=n ) where n is the desired ",
				"number of cores." )
				
			} else{
				
				options( mc.cores=n.cores )
				message( "\n'mc.cores' option set to ",
				as.character( n.cores ),".\n" )
				
			}
			
		}
		
		message( "Using ", as.character( n.cores ),
		" cores..." )
	
	}
			
	mean.shift.algorithm <- meanShiftAlgorithmAll( X=X, h=h,
	kernel=kernel, tol.stop=tol.stop, multi.core=multi.core, legacy.mode=legacy.mode)
	
	## find connected components
	message( "Finding clusters..." )
	output <- connectedComponents( X=mean.shift.algorithm,
	tol.epsilon=tol.epsilon )
	
	invisible( output )
	
}
