#include <RcppArmadillo.h>
// [[Rcpp::depends(Matrix,RcppArmadillo)]]
#include <iostream>
#include <armadillo>
#include <math.h>


using namespace Rcpp;
using namespace std;
using namespace arma;

//' C++ implementation of time consuming loops for plot
//'
//' @export
// [[Rcpp::export]]
arma::mat resCpp(Function fun, arma::vec theta, arma::vec s2)
{
  arma::vec w=as<arma::vec>(fun(theta,s2));
  double Dim = w.size();
  arma::mat res=randu<arma::mat>(Dim,100);
  for (int i=0; i<100; i++)
  {
    res.col(i) = as<arma::vec>(fun(theta,s2));
  }
  return res;
}



//' C++ implementation of time consuming loops for plot
//'
//' @export
// [[Rcpp::export]]
arma::mat resCppD(Function fun, arma::vec theta, arma::vec thetaD, arma::vec s2)
{
  arma::vec w=as<arma::vec>(fun(theta,thetaD,s2));
  double Dim = w.size();
  arma::mat res=randu<arma::mat>(Dim,100);
  for (int i=0; i<100; i++)
  {
    res.col(i) = as<arma::vec>(fun(theta,thetaD,s2));
  }
  return res;
}
