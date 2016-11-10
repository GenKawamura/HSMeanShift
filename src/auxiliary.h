#define KERNEL_1 "gaussianKernel"
#define KERNEL_FUNC_1 RcppGaussianKernel

#define KERNEL_2 "exponentialKernel"
#define KERNEL_FUNC_2 RcppExponentialKernel

#define KERNEL_3 "cubicKernel"
#define KERNEL_FUNC_3 RcppCubicKernel

#define KERNEL_4 "epanechnikovKernel"
#define KERNEL_FUNC_4 RcppEpanechnikovKernel

Rcpp::NumericVector RcppGaussianKernel(Rcpp::NumericVector x);
Rcpp::NumericVector RcppExponentialKernel(Rcpp::NumericVector x);
Rcpp::NumericVector RcppCubicKernel(Rcpp::NumericVector x);
Rcpp::NumericVector RcppEpanechnikovKernel(Rcpp::NumericVector x);
double RcppDistanceFunction(Rcpp::NumericVector x, Rcpp::NumericVector y );


