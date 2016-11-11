findClosestLargerPowerOf2 <- function( x ){
	
	if( x > 0 ){
		
		out <- ceiling( log( x, 2 ) )
		return( out )
		
	} else{
		
		stop( "x must be a positive number." )
		invisible( NULL )
		
	}
	
}

###

normalizeToUnitInterval <- function( x ){
	
	## function to normalize the domain of an observed
	## curve to the unit interval
	
	range.x <- range( x )
	
	output <- ( x - range.x[1] ) / ( range.x[2] - range.x[1] )
	
	return( output )
	
}

#' Function to project a curve on a wavelet basis.
#'
#' This function performs the Discrete Wavelet Transform (DWT) on a numeric vector representing a curve (i.e. a "functional" datum) observed on a grid and thresholds the wavelet coefficients, thus yielding a denoised and compressed representation of the same curve.
#'
#' The function normalizes the input grid to the standard unit interval, i.e. the minimum and the maximum values of \code{x.grid} are respectively 0 and 1.
#'
#' \code{projectCurveWavelet} is designed to be used as a preliminary step towards functional clustering using the mean shift algorithm. Given a sample of curves, \code{projectCurveWavelet} can be used to represent each curve as a sparse vector of coefficients. These coefficients can be fed as a matrix to \code{\link{msClustering}} or \code{\link{bmsClustering}} and clustered via the mean shift algorithm or the blurring mean shift algorithm.
#'
#'
#' @param x a numeric vector of x coordinates at which the curve is observed.
#' @param y a numeric vector of y coordinates representing the curve. \code{x} and \code{y} must have the same length.
#' @param irreg.grid logical. TRUE if \code{x} is not an equispaced grid.
#' @param grid.length a positive power of 2 or NULL (default). In order to apply the DWT, \code{length(x)} must be a positive power of 2. By default, if \code{grid.length=NULL} and \code{length(x)} is not a power of 2, \code{x} is extended to an equispaced grid whose length is positive power of 2 and \code{y} is extended interpolated on the extended grid. If \code{projectCurveWavelets} is used on multiple curves, \code{grid.length} should be set manually to ensure that all the discretized curves have the same length before the DWT is applied on each of them.
#' @param filter.number an integer specifying the smoothness of the wavelet used in the wavelet decomposition of \code{y}. See the functions \code{\link{wd}} and \code{\link{irregwd}} of \pkg{wavethresh} for details.
#' @param family a character string specifying the family of wavelets used in the wavelet decomposition of \code{y}. See the functions \code{\link{wd}} and \code{\link{irregwd}} of \pkg{wavethresh} for details.
#' @param bc a character string specifying how to handle the boundary condition. See the functions \code{\link{wd}} and \code{\link{irregwd}} of \pkg{wavethresh} for details.
#' @param verbose logical. Controls the printing of "informative" messages whilst the computation progresses. Such messages are generally annoying so it is turned off by default.
#' @param ... further arguments to control the thresholding of the wavelet coefficients. See \code{\link{threshold.wd}} and \code{\link{threshold.irregwd}} of the \pkg{wavethresh} package for details. By default, \code{projectCurvesWavelets} uses the default values of \code{\link{threshold.wd}} and \code{\link{threshold.irregwd}} to perform the thresholding of the wavelet coefficients.
#'
#'
#' @return The function outputs a list with names
#'  \item{coefficients }{a numeric vector of thresholded wavelet coefficients.}
#'  \item{y.wdT }{an object of class \code{wd} or \code{irregwd}. See \code{\link{threshold.wd}} and \code{\link{threshold.irregwd}} of the \pkg{wavethresh} package for details.}
#'  \item{y.wavelet }{a numeric vector with the reconstruction of \code{y} after the application of the DWT and the thresholding of the wavelet coefficients.}
#'  \item{x.grid }{the extended and equispaced grid of x values associated to \code{y.wavelet}.}
#'
#'
#' @example examples/example-projectCurveWavelets.R
#'
#' @author Mattia Ciollaro and Daren Wang
#' @references Nason, G. (2010) \emph{Wavelet methods in statistics with R}.
#' @seealso \code{\link{wavethresh}} \code{\link{wd}} \code{\link{irregwd}} \code{\link{threshold.wd}} \code{\link{threshold.irregwd}} \code{\link{wr}} \code{\link{msClustering}} \code{\link{bmsClustering}}
#'
#' @useDynLib MeanShift
#' @export
#'
projectCurveWavelets <- function( x, y, irreg.grid=FALSE, grid.length=NULL,
filter.number=10, family="DaubLeAsymm", bc="periodic", verbose=FALSE, ... ){
	
		## REQUIRES 'wavethresh' PACKAGE
	
		## normalize to unit interval
		x <- normalizeToUnitInterval( x )
		length.x <- length( x )
		
		## irregular grid?
		if( irreg.grid ){
			
			if( is.null( grid.length ) ){
				
				if( log( length.x, 2 ) %% 1 == 0 ){
					
					grid.length <- length.x
					
				} else{
					
					closest.power <- findClosestLargerPowerOf2( length.x )
					grid.length <- 2^closest.power
								
				}
					
			} else if( grid.length <= 0 || ( log( grid.length, 2 ) %% 1 != 0 ) ){

				stop( "grid.length must be a positive power of 2.")				
				
			}
				
			## make regular grid
			grid <- makegrid( x, y, gridn=grid.length )
			x <- grid$gridt
			y <- grid$gridy
				
			## wavelet transform
			y.wd <- irregwd( gd=grid, filter.number=filter.number, family=family, 
			bc=bc, verbose=verbose )					
				
		} else if( is.null( grid.length ) ){
			
			if( log( length.x, 2 ) %% 1 == 0  ){
				
				## wavelet transform
				y.wd <- wd( y, filter.number=filter.number, family=family, 
				bc=bc, verbose=verbose )
				
			} else{
				
				closest.power <- findClosestLargerPowerOf2( length.x )
				grid.length <- 2^closest.power

				## make regular grid
				grid <- makegrid( x, y, gridn=grid.length )
				x <- grid$gridt
				y <- grid$gridy

				## wavelet transform
				y.wd <- irregwd( gd=grid, filter.number=filter.number, family=family, 
				bc=bc, verbose=verbose )					
				
			}
		
		} else if( grid.length <= 0 || ( log( grid.length, 2 ) %% 1 != 0 ) ){
			
				stop( "grid.length must be a positive power of 2.")							
			
		} else{
			
			## make regular grid
			grid <- makegrid( x, y, gridn=grid.length )
			x <- grid$gridt
			y <- grid$gridy			
			
			## wavelet transform
			y.wd <- irregwd( gd=grid, filter.number=filter.number, family=family, 
			bc=bc, verbose=verbose )					
			
		}
						
		## thresholding
		y.wdT <- threshold( y.wd, ... )
		
		## extract coefficients
		y.coeffs <- y.wdT$D
		
		## smooth curve
		y.wave <- wr( y.wdT )
		
		output <- list( coefficients=y.coeffs, y.wdT=y.wdT, y.wavelet=y.wave,
		x.grid=x )
		
		return( output )
		
}
