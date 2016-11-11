\name{bmsClustering}
\alias{bmsClustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform clustering using the blurring version of the mean shift algorithm.
}
\description{
This function implements the blurring mean shift algorithm, which approximates
the standard mean shift algorithm. Because it recursively updates the entire
sample at each iteration, the blurring version of the mean shift algorithm is
often faster than the standard version (especially if the standard mean shift algorithm is run using a single core).
}
\usage{
bmsClustering(X, h = NULL, kernel = "epanechnikovKernel",
tol.stop = 1e-06, max.iter = 100, tol.epsilon = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a \eqn{p \times n} matrix containing \eqn{n \ge 1} \eqn{p}-dimensional numeric vectors stored as columns. Each column of \code{X} represents a sample point.
}
  \item{h}{
a strictly positive bandwidth parameter.
}
  \item{kernel}{
a kernel function (as a character string). The following kernels are supported: 
\itemize{
\item Epanechnikov: \eqn{ K(x) = \frac{3}{2}(1-x^2)I_{[0,1]}(x) }; \code{kernel="epanechnikovKernel"}
\item cubic: \eqn{ K(x) = 4(1-x)^3I_{[0,1]}(x) }; \code{kernel="cubicKernel"}
\item Gaussian: \eqn{ K(x) = \sqrt{\frac{2}{\pi}}e^{-\frac{x^2}{2}}I_{[0,\infty)}(x) }; \code{kernel="gaussianKernel"}
\item exponential \eqn{ K(x) = e^{-x}I_{[0,\infty)}(x) }; \code{kernel="exponentialKernel"}.
}
The use of the Epanechnikov kernel is recommended when using the blurring
version of the mean shift algorithm.
}
  \item{tol.stop}{
a strictly positive tolerance parameter. The mean shift algorithm stops when its update generates a step of length smaller than \code{tol.stop}. \code{tol.stop} should be considerably smaller than \code{tol.epsilon}.
}
  \item{max.iter}{
a strictly positive integer specifying the maximum number of iterations before the algorithm is forced to stop.
}
  \item{tol.epsilon}{
a strictly positive tolerance parameter. Points that are less than \code{tol.epsilon}- separated are grouped in the same cluster once the algorithm stops.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
It is generally recommended to standardize \code{X} so that each variable has
unit variance prior to running the algorithm on the data.

Roughly speaking, larger values of \code{h} produce a coarser clustering (i.e. few and large clusters). For sufficiently large values of \code{h}, the algorithm produces a unique cluster containing all the data points. Smaller values of \code{h} produce a finer clustering (i.e. many small clusters). For sufficiently small values of \code{h}, each cluster that is identified by the algorithm will contain exactly one data point.

If \code{h} is not specified in the function call, then \code{h} is by default set to the 30th percentile of the empirical distribution of distances between the columns of \code{X}, i.e. \code{h=quantile( dist( t( X ) ), 0.3 )}.

In their implementation, \code{gaussianKernel} and \code{exponentialKernel} are rescaled to assign probability of at least 0.99 to the unit interval \eqn{[0,1]}. This ensures that all the kernels are roughly on the same scale.

When using the blurring version of the mean shift algorithm, it is generally recommended to use a compactly supported kernel. In particular, the algorithm is guaranteed to converge in finitely many iterations with the Epanechnikov kernel.
}
\value{
The function invisibly returns a list with names
\item{components }{a matrix containing the modes/cluster representatives by column.}
\item{labels }{an integer vector of cluster labels.}
}
\references{
Carreira-Perpinan, M. A. (2015) \emph{A review of mean-shift algorithms for clustering}. 
arXiv \url{http://arxiv.org/abs/1503.00687}
}
\author{
Mattia Ciollaro and Daren Wang.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{msClustering}}
}
\examples{
## an example using the iris dataset
## help( iris )

## prepare data matrix
iris.data <- t( iris[,c( "Sepal.Length", "Sepal.Width" )] )

## run blurring mean shift algorithm
clustering <- bmsClustering( iris.data )
print( clustering )

## plot the clusters
\dontrun{
plot( iris.data[1,], iris.data[2,], col=clustering$labels+2, cex=0.8,
pch=16, xlab="Sepal.Length", ylab="Sepal.Width" )
points( clustering$components[1,], clustering$components[2,],
col=2+( 1:ncol( clustering$components ) ), cex=1.8, pch=16 )}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
