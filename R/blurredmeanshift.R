blurringMeanShiftAlgorithm <- function( X, h=NULL,
kernel="epanechnikovKernel", tol.stop=1e-6, max.iter=100, legacy.mode=FALSE){
	
	if( is.null( h ) ){
		
		h <- quantile( dist( t( X ) ), 0.3 )
		
	}

        not.converged <- FALSE

        if ( legacy.mode ){

          close.enough <- FALSE

          old.X <- X

          iter <- 0
          
          ## while the largest update corresponds to a shift
          ## larger than 'tol.stop' and while number of iterations
          ## is smaller than 'max.iter'
          while( !close.enough ){
            
            ## apply blurring mean-shift operator and update
            iter <- iter + 1
            
            new.X <- blurringMeanShiftOperator( X=old.X, h=h, kernel=kernel )
            
            distance <- max( sqrt( colSums( old.X - new.X )^2 ) )
            
            old.X <- new.X
            
            close.enough <- ( distance < tol.stop )	
            
            if( iter >= max.iter ){
              
              not.converged <- TRUE
              break
              
            }
          }

	} else {
          ## Rcpp implementation of blurringMeanShiftAlgorithm. This method is much faster than ordinary R code.
          ## output <- RcppBlurringMeanShiftAlgorithm( X, h, kernel, tol.stop, max.iter )
          not.converged <- output[["not_converged"]]
          new.X <- output[["new_x"]]
        }
        
	
	if( not.converged ){
		
		if( kernel == "epanechnikovKernel"){
			
			warning( "Reached maximum number of iterations (", 
			as.character( max.iter),"). The algorithm ",
			"didn't converge. Try increasing max.iter." )
			
		} else{

			warning( "Reached maximum number of iterations (", 
			as.character( max.iter),"). The algorithm ",
			"didn't converge. Try kernel=\"epanechnikovKernel\"." )
			
		}
		
	} else {

		message( "Blurring mean-shift algorithm ran successfully.\n")
			
	}
	
	return( new.X )
	
}


#' Function to perform clustering using the blurring version of the mean shift algorithm.
#'
#' This function implements the blurring mean shift algorithm, which approximates the standard mean shift algorithm. Because it recursively updates the entire sample at each iteration, the blurring version of the mean shift algorithm is often faster than the standard version (especially if the standard mean shift algorithm is run using a single core).
#'
#' It is generally recommended to standardize \code{X} so that each variable has unit variance prior to running the algorithm on the data.
#'
#' Roughly speaking, larger values of \code{h} produce a coarser clustering (i.e. few and large clusters). For sufficiently large values of \code{h}, the algorithm produces a unique cluster containing all the data points. Smaller values of \code{h} produce a finer clustering (i.e. many small clusters). For sufficiently small values of \code{h}, each cluster that is identified by the algorithm will contain exactly one data point.
#'
#' If \code{h} is not specified in the function call, then \code{h} is by default set to the 30th percentile of the empirical distribution of distances between the columns of \code{X}, i.e. \code{h=quantile( dist( t( X ) ), 0.3 )}.
#' 
#' In their implementation, \code{gaussianKernel} and \code{exponentialKernel} are rescaled to assign probability of at least 0.99 to the unit interval \eqn{[0,1]}. This ensures that all the kernels are roughly on the same scale.
#'
#' When using the blurring version of the mean shift algorithm, it is generally recommended to use a compactly supported kernel. In particular, the algorithm is guaranteed to converge in finitely many iterations with the Epanechnikov kernel.
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
#' The use of the Epanechnikov kernel is recommended when using the blurring version of the mean shift algorithm.
#' @param tol.stop a strictly positive tolerance parameter. The mean shift algorithm stops when its update generates a step of length smaller than \code{tol.stop}. \code{tol.stop} should be considerably smaller than \code{tol.epsilon}.
#' @param max.iter a strictly positive integer specifying the maximum number of iterations before the algorithm is forced to stop.
#' @param tol.epsilon a strictly positive tolerance parameter. Points that are less than \code{tol.epsilon}- separated are grouped in the same cluster once the algorithm stops.
#' @param legacy.mode default is \code{FALSE}. If \code{TRUE}, legacy kernel codes written in R are used. The legacy R implementation is much slower.
#'
#'
#' @return The function invisibly returns a list with names
#' \item{components}{a matrix containing the modes/cluster representatives by column.}
#' \item{labels}{an integer vector of cluster labels.}
#'
#'
#' @example examples/example-bmsClustering.R
#'
#'
#' @author Mattia Ciollaro, Daren Wang and Gen Kawamura
#' @references Carreira-Perpinan, M. A. (2015) \emph{A review of mean-shift algorithms for clustering}. arXiv \url{http://arxiv.org/abs/1503.00687}
#' @seealso \code{\link{msClustering}}
#'
#' @useDynLib MeanShift
#' @export
#'
bmsClustering <- function( X, h=NULL, kernel="epanechnikovKernel",
tol.stop=1e-6, max.iter=100, tol.epsilon=1e-3, legacy.mode=FALSE ){
	
	# minimal input checking
	X <- as.matrix( X )
	max.iter <- as.integer( max.iter )
	
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
	
	if( max.iter <= 0 ){
		
		stop( "The maximum number of iterations must be a positive ",
		"integer." )
		
	}
	
	if( tol.stop <= 0 || tol.epsilon <= 0 ){
		
		stop( "All tolerances must be strictly positive.")
		
	}
	
	## run blurring mean-shift algorithm
	message( "\nRunning blurring mean-shift algorithm...\n" )
	
	blurring.mean.shift.algorithm <- blurringMeanShiftAlgorithm( X=X,
	h=h, kernel=kernel, tol.stop=tol.stop, max.iter=max.iter, legacy.mode=legacy.mode)
	
	## find connected components
	message( "Finding clusters..." )
	output <- connectedComponents( X=blurring.mean.shift.algorithm,
	tol.epsilon=tol.epsilon )
	
	invisible( output )

}
