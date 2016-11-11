\name{projectCurveWavelets}
\alias{projectCurveWavelets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to project a curve on a wavelet basis.
}
\description{
This function performs the Discrete Wavelet Transform (DWT) on a numeric vector representing a curve (i.e. a "functional" datum) observed on a grid and thresholds the wavelet coefficients, thus yielding a denoised and compressed representation of the same curve.
}
\usage{
projectCurveWavelets( x, y, irreg.grid=FALSE, grid.length=NULL,
filter.number=10, family="DaubLeAsymm", bc="periodic", verbose=FALSE, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector of x coordinates at which the curve is observed.
}
  \item{y}{
a numeric vector of y coordinates representing the curve. \code{x} and \code{y} must have the same length.
}
  \item{irreg.grid}{
  logical. TRUE if \code{x} is not an equispaced grid.
}
  \item{grid.length}{
  a positive power of 2 or NULL (default). In order to apply the DWT, \code{length(x)} must be a positive power of 2. By default, if \code{grid.length=NULL} and \code{length(x)} is not a power of 2, \code{x} is extended to an equispaced grid whose length is positive power of 2 and \code{y} is extended interpolated on the extended grid. If \code{projectCurveWavelets} is used on multiple curves, \code{grid.length} should be set manually to ensure that all the discretized curves have the same length before the DWT is applied on each of them.
}
  \item{filter.number}{
an integer specifying the smoothness of the wavelet used in the wavelet decomposition of \code{y}. See the functions \code{\link{wd}} and \code{\link{irregwd}} of \pkg{wavethresh} for details.
}
  \item{family}{
a character string specifying the family of wavelets used in the wavelet decomposition of \code{y}. See the functions \code{\link{wd}} and \code{\link{irregwd}} of \pkg{wavethresh} for details.
}
  \item{bc}{
a character string specifying how to handle the boundary condition. See the functions \code{\link{wd}} and \code{\link{irregwd}} of \pkg{wavethresh} for details.
}
  \item{verbose}{
logical. Controls the printing of "informative" messages whilst the computation progresses. Such messages are generally annoying so it is turned off by default.
}
  \item{...}{
further arguments to control the thresholding of the wavelet coefficients. See \code{\link{threshold.wd}} and \code{\link{threshold.irregwd}} of the \pkg{wavethresh} package for details. By default, \code{projectCurvesWavelets} uses the default values of \code{\link{threshold.wd}} and \code{\link{threshold.irregwd}} to perform the thresholding of the wavelet coefficients.
}
}
\details{
The function normalizes the input grid to the standard unit interval, i.e. the minimum and the maximum values of \code{x.grid} are respectively 0 and 1.

\code{projectCurveWavelet} is designed to be used as a preliminary step towards functional clustering using the mean shift algorithm. Given a sample of curves, \code{projectCurveWavelet} can be used to represent each curve as a sparse vector of coefficients. These coefficients can be fed as a matrix to \code{\link{msClustering}} or \code{\link{bmsClustering}} and clustered via the mean shift algorithm or the blurring mean shift algorithm.
}
\value{
The function outputs a list with names
\item{coefficients }{a numeric vector of thresholded wavelet coefficients.}
\item{y.wdT }{an object of class \code{wd} or \code{irregwd}. See \code{\link{threshold.wd}} and \code{\link{threshold.irregwd}} of the \pkg{wavethresh} package for details.}
\item{y.wavelet }{a numeric vector with the reconstruction of \code{y} after the application of the DWT and the thresholding of the wavelet coefficients.}
\item{x.grid }{the extended and equispaced grid of x values associated to \code{y.wavelet}.}
}
\references{
Nason, G. (2010) \emph{Wavelet methods in statistics with R}.
}
\author{
Mattia Ciollaro and Daren Wang
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{wavethresh}} \code{\link{wd}} \code{\link{irregwd}} \code{\link{threshold.wd}} \code{\link{threshold.irregwd}} \code{\link{wr}} \code{\link{msClustering}} \code{\link{bmsClustering}}
}
\examples{
## generate a noisy curve observed on a regular grid
set.seed( 1 )
n.grid <- 1000
x <- seq( 2, 8, length=n.grid )
sigma.epsilon1 <- 2
sigma.epsilon2 <- 2.5
sigma.epsilon3 <- 3
sigma.epsilon4 <- 1
epsilon <- rnorm( 1000, sd=rep( c( sigma.epsilon1,
sigma.epsilon2, sigma.epsilon3, sigma.epsilon4 ),
rep( 250, 4 ) ) )
y <- x*sin( 3*x ) + 0.3*x^2 + epsilon

## project on wavelet basis with soft universal thresholding
## of the wavelet coefficients
wave <- projectCurveWavelets( x, y, type="soft", policy="universal" )

## plot wavelet reconstruction of the curve
\dontrun{
x.norm <- ( x - min( x ) ) / ( max( x ) - min( x ) )
plot( x.norm, y )
lines( wave$x.grid, wave$y.wavelet, col=2, lwd=3 )}

## inspect wavelet coefficients
wave.coeffs <- wave$coefficients
print( length( wave.coeffs ) )  ## 1023 coefficients
print( sum( wave.coeffs != 0 ) )  ## only 12 are non-zero
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
