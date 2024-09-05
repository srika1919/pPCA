## Compute a Linear Combination of Matrix-Vector Products
##
## @author Dr.Somak Dutta, Srika Raja
## @decription
## This function calculates a linear combination of matrix-vector products based on
## the provided input vector and a list of matrices.
##
## @param v A numeric vector to be multiplied with the matrices.
## @param args A list containing the following elements:
##   \item{X}{A list of matrices to be multiplied with segments of v.}
##   \item{m}{An integer specifying the number of matrices in X.}
##   \item{cns}{A numeric vector of cumulative sums, used to determine
##              the segments of v to be used with each matrix.}
##
## @return A numeric vector resulting from the sum of matrix-vector products.



Afun = function(v, args) {
  X = args$X
  m = length(X)
  cns = args$cns
  u = 0
  for(j in 1:m)
    u = u + X[[j]] %*% v[ {cns[j]+1}:cns[j+1] ]

  if(class(u)[1] == "dgeMatrix") {
    return(u@x)
  } else if(is(u,"numeric")) {
    return(u)
  } else
    return(as(u,"numeric"))
}
