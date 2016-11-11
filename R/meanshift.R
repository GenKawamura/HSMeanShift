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
#'
#' This function implements the mean shift algorithm. The algorithm locates the modes of a kernel density estimator and associates each data point to exactly one of the modes, thus effectively clustering the data.
#'
#'
#' It is generally recommended to standardize \code{X} so that each variable has unit variance prior to running the algorithm on the data.
#'
#' Roughly speaking, larger values of \code{h} produce a coarser clustering (i.e. few and large clusters). For sufficiently large values of \code{h}, the algorithm produces a unique cluster containing all the data points. Smaller values of \code{h} produce a finer clustering (i.e. many small clusters). For sufficiently small values of \code{h}, each cluster that is identified by the algorithm will contain exactly one data point.
#'
#' If \code{h} is not specified in the function call, then \code{h} is by default set to the 30th percentile of the empirical distribution of distances between the columns of \code{X}, i.e. \code{h=quantile( dist( t( X ) ), 0.3 )}.
#'
#' In their implementation, \code{gaussianKernel} and \code{exponentialKernel} are rescaled to assign probability of at least 0.99 to the unit interval \eqn{[0,1]}. This ensures that all the kernels are roughly on the same scale.
#' 
#' To specify the number of cores when \code{multi.core=TRUE}, the option
#' \code{mc.cores} needs to be set with \code{options( mc.cores=n.cores )}, where \code{n.cores} is the number of cores that the mean shift algorithm is allowed to use for parallel computation.
#'
#'
#' @param X a \eqn{p \times n} matrix containing \eqn{n \ge 1} \eqn{p}-dimensional numeric vectors stored as columns. Each column of \code{X} represents a sample point.
#' @param h a strictly positive bandwidth parameter.
#' @param kernel a kernel function (as a character string). The following kernels are supported:
#' \itemize{
#'  \item Epanechnikov: \eqn{ K(x) = \frac{3}{2}(1-x^2)I <- {[0,1]}(x) }; \code{kernel="epanechnikovKernel"}
#'  \item cubic: \eqn{ K(x) = 4(1-x)^3I <- {[0,1]}(x) }; \code{kernel="cubicKernel"}
#'  \item Gaussian: \eqn{ K(x) = \sqrt{\frac{2}{\pi}}e^{-\frac{x^2}{2}}I <- {[0,\infty)}(x) }; \code{kernel="gaussianKernel"}
#'  \item exponential \eqn{ K(x) = e^{-x}I <- {[0,\infty)}(x) }; \code{kernel="exponentialKernel"}.
#' }
#' @param tol.stop a strictly positive tolerance parameter. The algorithm stops when all of the updates generate steps of length smaller than \code{tol.stop}. \code{tol.stop} should be considerably smaller than \code{tol.epsilon}.
#' @param tol.epsilon a strictly positive tolerance parameter. Points that are less than \code{tol.epsilon}- separated are grouped in the same cluster once the algorithm stops.
#' @param multi.core logical. If \code{TRUE}, the mean shift algorithm is parallelized.
#' @param legacy.mode default is \code{FALSE}. If \code{TRUE}, legacy kernel codes written in R are used. The legacy R implementation is much slower.
#'
#'
#' @return The function invisibly returns a list with names
#' \item{components}{a matrix containing the modes/cluster representatives by column.}
#' \item{labels}{an integer vector of cluster labels.}
#'
#'
#' @example examples/example-msClustering.R
#'
#' @author Mattia Ciollaro, Daren Wang and Gen Kawamura
#' @references Carreira-Perpinan, M. A. (2015) \emph{A review of mean-shift algorithms for clustering}. arXiv \url{http://arxiv.org/abs/1503.00687}
#' @seealso \code{\link{bmsClustering}}
#'
#' @useDynLib MeanShift
#' @export
#'
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
