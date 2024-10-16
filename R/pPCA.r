#' Performs a principal component analysis on a large sparse matrices or a list of large sparse
#' matrices and returns the results as an object compatible to class prcomp
#'
#' @author Srika Raja and Somak Dutta
#'
#' @description
#' Performs a partial principal component analysis on a large sparse matrices or a list of large sparse
#' matrices and returns the results as an object compatible to class prcomp. Uses RSpectra library
#' to compute the largest eigenvalues.
#'
#' @param x A matrix, sparse matrix (Matrix::dgCMatrix), or a list of these. When a list
#' is supplied, the entries are concatenated horizontally (implicitly). See description.
#' @param rank An integer specifying the number of principal components to compute.
#' @param retX A logical value indicating whether the rotated variables (PC scores) should be returned.
#' @param scale. A logical value indicating whether the variables should be scaled to have
#' unit variance before the analysis takes place.
#' @param normalize A logical value indicating whether the principal component scores should be normalized.
#' @param sd.tol A positive number, warnings are printed if the standard deviation of any
#' column is less than this threshold.
#'
#' @return pPCA returns a list with class "pPCA" (compatible with "prcomp") containing the following
#' components:
#'   \item{sdev}{A vector of the singular values (standard deviations of the principal components).}
#'   \item{rotation}{A matrix whose columns contain the eigenvectors (loadings).}
#'   \item{x}{A matrix of the principal component scores, returned if retX is true. This is
#'   the centred (and scaled if requested) data multiplied by the rotation matrix.}
#'   \item{center}{column means.}
#'   \item{scale}{column standard deviations, if scale. is true. Otherwise, FALSE.}
#'
#' @details
#' When the input argument is a matrix (of class "matrix" or "dgCMatrix"), principal component analysis
#'  is performed to extract a few largest components. When a list of matrices is passed, the partial PCA
#'  is performed on the horizontally concatenated matrix, i.e., if \code{x = list(X1,X2,X3)} then the
#'  partial PCA is done on the matrix [X1 X2 X3], without concatenating the matrices explicitly. This can be
#'  useful when the matrix is so high-dimensional that the total number of non-zero entries
#'  exceed 2^31-1 (roughly 9.33e10), the capacity of a 32 bit integer. For example, in PCA with very
#'  high-dimensional SNP data, the sparse matrices can be stored for each chromosome within the capacity
#'   of 32 bit integers.
#'
#' @note
#' The partial SVD is computed through the RSpectra package. All elements in the first row of the rotation
#' matrix are positive.
#' @references Raja, S. and Dutta, S. (2024). Matrix-free partial PCA of partitioned genetic data.
#' REU project 2024, Iowa State University.
#' @references Dai, F., Dutta, S., and, Maitra, R. (2020). A Matrix-Free Likelihood Method for
#'Exploratory Factor Analysis of High-Dimensional Gaussian Data. Journal of Computational and
#'Graphical Statistics, 29(3), 675--680.
#' @seealso \code{\link{biplot},\link{prcomp}}




#'
#' @examples
#'
#' library(Matrix)
#' set.seed(20190329)
#' m <- rsparsematrix(50,100,density = 0.35)
#' results <- pPCA(m, rank = 2)
#' biplot(results)
#' data <- list(rsparsematrix(nrow = 50,ncol = 10,density = 0.35),
#'              rsparsematrix(nrow = 50,ncol = 40,density = 0.35)) # Using a list of matrices
#' result <- pPCA(data, rank = 3)
#' print(result)
#' biplot(result)
#'
#'

pPCA <- function(x, rank, retX = TRUE, scale. = TRUE, normalize = FALSE, sd.tol = 1e-5) {

  if(is.list(x))
    result <- pca_list(x,rank,retX,scale., normalize, sd.tol)
  else
    result <- pca_matrix(x,rank,retX,scale., normalize, sd.tol)

  return(result)
}
