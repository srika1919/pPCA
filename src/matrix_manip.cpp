/*
 * Author: Somak Dutta
 */

#include <Rcpp.h>


using namespace Rcpp;



// [[Rcpp::export]]
NumericVector colMSD_dgc(S4 mat,NumericVector m) {
  IntegerVector dims = mat.slot("Dim");
  int ncol = dims[1];
  int nrow = dims[0];

  int *p = INTEGER(mat.slot("p"));
  double *x = REAL(mat.slot("x"));

  NumericVector y(ncol);
  for(int j=0;j<ncol;++j)
  {
    double ssq = 0.0;
    double mj = m[j];
    double *start = &x[p[j]];
    double *end = &x[p[j+1]];
    for(;start < end; ++start)
    {
      double temp = *start - mj;
      ssq += temp*temp;
    }
    y[j] = (ssq + mj*mj*(nrow - (p[j+1] - p[j] + 0.0)))/(nrow-1);
  }
  return(y);
}


