#include <Rcpp.h>
#include <R.h>
#include <math.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
NumericVector tnorm(int n, double lo, double hi,
                    NumericMatrix mu, NumericMatrix sig) {
  int musize = mu.nrow()*mu.ncol();
  NumericVector lotemp;
  if (musize > 1) {
    NumericVector temp(lo, musize);
    lotemp = temp;
  }
  NumericVector hitemp;
  if (musize > 1) {
    NumericVector temp(hi, musize);
    hitemp = temp;
  }
  
  NumericVector q1 = R::pnorm(lotemp, mu, sig, true, false);
  NumericVector q2 = R::pnorm(hitemp, mu, sig, true, false);
  
  NumericVector z = R::runif(n, q1, q2);
  z = R::qnorm(z, mu, sig);
  return z;
}

// [[Rcpp::export]]
NumericVector fixMH(NumericVector lo, NumericVector hi,
                    NumericVector old1, NumericVector new1, NumericVector jump) {
  NumericVector jold = R::pnorm(hi, old1, jump, true, false) -
    R::pnorm(lo, old1, jump, true, false);
  NumericVector jnew = R::pnorm(hi, new1, jump, true, false) -
    R::pnorm(lo, new1, jump, true, false);
    return log(jold) - log(jnew);//vetor?????
}

// [[Rcpp::export]]
NumericVector getLogl(NumericVector theta, double phi, NumericVector y, int nm) {
  NumericVector prob(theta.length());
  for (int i=0; i< theta.length(); i++)
    prob[i] = theta[i]*phi;
  for (int i=0; i<prob.length(); i++)
    if (prob[i] < 1e-05)
      prob[i] = 1e-05;
  
  NumericVector cond(prob.length());
  for (int i=0; i<prob.length(); i++)
    if(prob[i]>0.99999)
      prob[i] = 0.99999;
  return R::dbinom(y, nm, prob, true);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
