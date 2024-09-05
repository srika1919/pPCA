## Compute Transposed Linear Combination of Matrix-Vector Products
##
## @author Dr.Somak Dutta, Srika Raja
##
## @description
## This function calculates a transposed linear combination of matrix-vector products based on the provided input vector and a list of matrices.
##
## @param v A numeric vector to be multiplied with the transposed matrices.
## @param args A list containing the following elements:
##   \item{X}{A list of matrices to be transposed and multiplied with v.}
##   \item{m}{An integer specifying the number of matrices in X.}
##   \item{cns}{A numeric vector of cumulative sums, used to determine
##              the segments of the output vector to be filled.}
##   \item{p}{An integer specifying the length of the output vector.}
##
##@return A numeric vector of length p, resulting from the transposed matrix-vector products


Atrans = function(v,args) {
  X = args$X
  m = length(X)
  cns = args$cns
  p = args$p
  u = numeric(p)
  for(j in 1:m)
    u[{cns[j]+1}:cns[j+1] ] = as.numeric(crossprod(X[[j]] , v))

  return(u)
}
