#include <Rcpp.h>
#include <cmath>
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
NumericMatrix GaussianKernel(NumericMatrix X, double sigma)
{
  int nrow = X.nrow();
  int ncol = X.ncol();
  double temp = 0;
  NumericMatrix result = NumericMatrix(nrow, nrow);
  for(int i=0; i< nrow; i++)
  {
    for(int j=0; j<nrow; j++)
    {
      temp = 0;
      for(int k=0; k<ncol; k++)
      {
        temp = temp + (X(i,k) - X(j,k)) * (X(i,k) - X(j,k));
      }
      temp = -temp / (sigma * sigma);
      result(i, j) = exp(temp);
    }
  }
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


