#include <RcppArmadillo.h>
#include "auxiliary.h"

using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix RcppBlurringMeanShiftOperator(const NumericMatrix & x, double h, String kernel){
  int n = x.ncol();
        
  // mean-shift operator
  // compute distances
  NumericVector distances(n);
  for (int i = 0; i < n; i++) distances[i] = RcppDistanceFunction(x(_,i), x);

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
    output = wrap(as<arma::mat>(x) * as<arma::vec>(kernel_weights));
  } else{
    NumericMatrix a = x;
  }
  NumericMatrix k = clone(x);
  return k;

}


// [[Rcpp::export]]
List RcppBlurringMeanShiftAlgorithm(const NumericMatrix & x, double h, String kernel, double tol_stop, int max_iter){

  bool not_converged = false;

  NumericMatrix new_x = clone(x);
  NumericMatrix old_x = clone(x);
  int iter = 0;
    
  // while the trajectory has not converged
  // (update produced a shift larger than tol_stop)
  while( true ) {
    // apply blurring mean-shift operator and update                                                                                                  
    iter++;

    new_x = RcppBlurringMeanShiftOperator(old_x, h, kernel);
    //double distance <- max( sqrt( colSums( old.X - new.X )^2 ) )
    double distance = 0;
    old_x = clone(new_x);

    // close enough
    if ( distance < tol_stop ) break;

    // reached number of max iteration
    if( iter >= max_iter ){
      not_converged = true;
      break;
    }
  }

  // output
  return List::create(Named("not_converged", not_converged), Named("new_x", new_x));
}

