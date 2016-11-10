#include <RcppArmadillo.h>
#include "auxiliary.h"

using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector RcppMeanShiftOperator(NumericVector x, const NumericMatrix & points, double h, String kernel){
  int n = points.ncol();
        
  // mean-shift operator
  // compute distances
  NumericVector distances(n);
  for (int i = 0; i < n; i++) distances[i] = RcppDistanceFunction(points(_,i), x);

  // scale by bandwidth
  NumericVector scaled_distances(n);
  for (int i = 0; i < n; i++) scaled_distances[i] = distances[i] / h;

  // evaluate kernel
  NumericVector kernel_values(n);
  if ( kernel == KERNEL_1 ) kernel_values = KERNEL_FUNC_1( scaled_distances );
  if ( kernel == KERNEL_2 ) kernel_values = KERNEL_FUNC_2( scaled_distances );
  if ( kernel == KERNEL_3 ) kernel_values = KERNEL_FUNC_3( scaled_distances );
  if ( kernel == KERNEL_4 ) kernel_values = KERNEL_FUNC_4( scaled_distances );

  // weights denominator
  double total_sum = sum(kernel_values);

  // initialization of variables
  NumericVector kernel_weights(n);
  NumericVector output(n);

  // mean-shift weights
  //Rprintf("total_sum = %e\n", total_sum);  
  if( total_sum > 0 ){
    // update
    for (int i = 0; i < n; i++) kernel_weights[i] = kernel_values[i] / total_sum;
    //output = points %*% kernel_weights;
    output = wrap(as<arma::mat>(points) * as<arma::vec>(kernel_weights));
  } else{
    output = x;
  }
  return output;

}


// [[Rcpp::export]]
NumericVector RcppMeanShiftAlgorithm(NumericVector x, const NumericMatrix & points, double h, String kernel, double tol_stop){

  NumericVector new_x(x.size());
  NumericVector old_x = clone(x);
    
  // while the trajectory has not converged
  // (update produced a shift larger than tol_stop)
  while( true ) {
    // apply mean-shift operator and update
    new_x = RcppMeanShiftOperator(old_x, points, h, kernel);
    double distance = RcppDistanceFunction( old_x, new_x );
    old_x = clone(new_x);

    // close enough
    if ( distance < tol_stop ) break;
  }

  return new_x;
}

