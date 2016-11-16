#include <math.h>
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;


double computeGaussianKernel(double y) {
  double value = 0;

  if (0 <= y){
    //value = (double) 2 / 0.388 * dnorm( y / 0.388 );
    NumericVector k(1);
    k = R::dnorm( y / 0.388, 0, 1, false);
    value = (double) 2 / 0.388 * k[0];
  }
  return value;
}


// [[Rcpp::export]]
NumericVector RcppGaussianKernel(NumericVector x){
  int n = x.size();
  NumericVector output(n);
  for (int i = 0; i < n; i++) output[i] = computeGaussianKernel(x[i]);

  return output;
}


double computeExponentialKernel(double y) {
  double value = 0;
  if (0 <= y){
    // value <- dexp( y, rate=4.61 )
    value = R::dexp(y, 4.61, false);
  }
  return value;
}


// [[Rcpp::export]]
NumericVector RcppExponentialKernel(NumericVector x){
  int n = x.size();
  NumericVector output(n);
  for (int i = 0; i < n; i++) output[i] = computeExponentialKernel(x[i]);

  return output;
}


double computeCubicKernel(double y) {
  double value = 0;
  if( 0 <= y && y<= 1 ){
    value = 4 * pow(1 - y, 3);
  }
  return value;
}


// [[Rcpp::export]]
NumericVector RcppCubicKernel(NumericVector x){
  int n = x.size();
  NumericVector output(n);
  for (int i = 0; i < n; i++) output[i] = computeCubicKernel(x[i]);

  return output;
}


double computeEpanechnikovKernel(double y) {
  double value = 0;
  if (0 <= y && y <= 1){
    value = (double) 3 / 2 * (1 - pow(y, 2));
  }
  return value;
}


// [[Rcpp::export]]
NumericVector RcppEpanechnikovKernel(NumericVector x){
  int n = x.size();
  NumericVector output(n);
  for (int i = 0; i < n; i++) output[i] = computeEpanechnikovKernel(x[i]);

  return output;
}


// [[Rcpp::export]]
double RcppDistanceFunction(NumericVector x, NumericVector y ){
	
  // function to compute the standard euclidean distance
  // output <- sqrt( sum( ( x - y )^2 ) )

  double output = 0;
  int n = x.size();
  // sum((x-y)^2)
  for (int i = 0; i < n; i++) output += pow(x[i] - y[i], 2);
  return sqrt(output);
}

