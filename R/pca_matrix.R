
## Perform Principal Component Analysis using Singular Value Decomposition for a sparse matrix.
##
## @author Srika Raja and Somak Dutta
##
## @description
## This function performs PCA on a matrices using singular value decomposition (SVD) and gives the top k values.
##
## @param m A sparse matrix to perform PCA.
## @param k An integer specifying the number of principal components to compute.
##
## @return A list containing three elements:
##   \item{sdev}{A vector of the singular values (standard deviations of the principal components).}
##   \item{rotation}{A matrix whose columns contain the eigenvectors (loadings).}
##   \item{pcs}{A matrix of the principal component scores.}
##
## @note
## This function assumes that the 'svds' function from the 'RSpectra⁠' package is available.
## Make sure to have the package installed and loaded before using this function.

## @examples
## # library(Matrix)
## # data <- rsparsematrix(25,30,density = 0.3)
## # result <- pca_sv(data, k = 3)
## # print(result)
##

pca_matrix <- function(m, rank, retX = TRUE, scale. = TRUE, normalize = FALSE, sd.tol = 1e-5){


  m_colnames <- colnames(m)
  m_rownames <- rownames(m)
  dimnames(m) <- list(NULL,NULL) # to avoid copying dimnaes unnecessarily during matrix-vector products
  # This does not create copies of the matrix thanks to read-on-write.

  cm <- colMeans(m)
  n <- nrow(m)
  p <- ncol(m)

  if(rank >= min(n,p))
    stop("Number of principal components cannot be more than any dimension")

  if(rank > min(n,p)/4)
    warning("Too many principal components requested.")


  if(class(m)[1] == "dgCMatrix") {
    sds <- sqrt(colMSD_dgc(m,cm))
  } else if(class(m)[1] == "matrix") {
    sds <- apply(m,2,sd)
  } else {
    stop("Only base::matrix and Matrix::dgCMatrix matrices (or a list of these) are supported.")
  }






  if(scale.) {
    if(any(sds == 0))
      stop("cannot rescale a constant/zero column to unit variance")
    if(any(sds < sd.tol))
      warning("Columns with very low sd (< ",sd.tol,") encountered. They should be removed",immediate. = T)
    sc = sqrt(n-1)*sds
  } else {
    sc = rep(sqrt(n-1),p)
  }



  sv <- svds(m, k=rank, nu=0, nv=rank, opts = list(center = cm,scale = sc))

  # Ensure first row of v is positive
  for(ii in 1:ncol(sv$v))
    if(sv$v[1,ii] < 0) sv$v[,ii] = -sv$v[,ii]

  dimnames(sv$v) <- list(m_colnames, paste0("PC", seq_len(ncol(sv$v))))

  if(retX)
  {
    if(scale.) {
      w <- sv$v / sds;
    } else {
      w = sv$v;
    }
    pcscores = sweep( as.matrix(m %*% w) ,MARGIN = 2,STATS =  t(cm) %*% w,FUN = "-")
    rownames(pcscores) <- m_rownames
  }
  result <- list("sdev" = sv$d, "rotation" =sv$v,center = cm)

  result$scale <- if(!scale.) FALSE else sds;
  result$sds <- sds;

  if(retX){
    if(normalize){
      result$x <- sweep(pcscores, 2, sv$d ,FUN = "/")
    }
    else{
      result$x <- pcscores
    }
  }

  result$nsample <- n

  class(result) <- c("pPCA","prcomp")
  return (result)
}



