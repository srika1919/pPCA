

## Perform Principal Component Analysis on a List of Matrices using Singular Value Decomposition
##
## @author Srika Raja and Somak Dutta
##
## @description
## This function performs PCA on a list of matrices using singular value decomposition (SVD).
##
## @param m A list of sparse matrices to perform PCA.
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
##
## @examples
## # library(Matrix)
## # data <- list(rsparsematrix(nrow = 25,ncol = 10,density = 0.35),rsparsematrix(nrow = 25,ncol = 40,density = 0.35))
## # result <- pca_sv_list(data, k = 3)
## # print(result)
##
##


pca_list <- function(lm, rank, retX = TRUE, scale. = TRUE, sd.tol = 1e-5){

  m_colnames <- lapply(lm,colnames)
  m_colnames <- if(any(sapply(m_colnames,is.null))) NULL else unlist(m_colnames)
  m_rownames <- rownames(lm[[1]])

  for(ii in 1:length(lm)) {
    dimnames(lm[[ii]]) <- list(NULL,NULL)
  }

  cm <- lapply(lm, colMeans)
  n <- sapply(lm,nrow)
  if(length(unique(n)) != 1) stop("Number of rows of each matrix in the list must be the same.")
  n <- n[1]

  ncols <- sapply(lm,ncol)
  cns <- c(0L,cumsum(ncols))
  p <- sum(ncols)

  if(rank >= min(n,p))
    stop("Number of principal components cannot be more than any dimension")

  if(rank > min(n,p)/4)
    warning("Too many principal components requested.")


  if(scale.) {
    sds = vector(mode = "list",length = length(lm))
    for(ii in 1:length(lm)) {
      m <- lm[[ii]]
      if(class(m)[1] == "dgCMatrix") {
        sds[[ii]] <- sqrt(colMSD_dgc(m,cm[[ii]]))
      } else if(class(m)[1] == "matrix") {
        sds[[ii]] <- apply(m,2,sd)
      } else {
        stop("Only base::matrix and Matrix::dgCMatrix matrices (or a list of these) are supported.")
      }
    }
    sds <- unlist(sds)
    if(any(sds == 0))
      stop("cannot rescale a constant/zero column to unit variance")

    if(any(sds < sd.tol))
      warning("Columns with very low sd (< ",sd.tol,") encountered. They should be removed",immediate. = T)

    sc = sqrt(n-1)*sds

  } else {
    sc = FALSE
  }

  cm <- unlist(cm)

  sv <- svds(A = Afun,k = rank,nu = 0,nv = rank,Atrans = Atrans,dim = c(n,p),
             args = list(X=lm,cns=cns,p=p),
             opts = list(center = cm,scale =sc))

  # Ensure first row of v is positive
  for(ii in 1:ncol(sv$v))
    if(sv$v[1,ii] < 0) sv$v[,ii] = -sv$v[,ii]

  dimnames(sv$v) <- list(m_colnames, paste0("PC", seq_len(ncol(sv$v))))

  # w <- matrix(0, ncol = rank, nrow = p)
  if(retX)
  {
    if(scale.) {
      w <- sv$v / sds;
    } else {
      w = sv$v;
    }
    pcscores <- matrix(0, nrow = n, ncol = rank)
    for (i in 1:ncol(sv$v)) {
      pcscores[, i] <- Afun(w[, i], args = list(X = lm, cns = cns, p = p))
    }
    pcscores <- sweep(pcscores, MARGIN =  2, STATS =  crossprod(cm,w),FUN = "-")
    rownames(pcscores) <- m_rownames
  }

  result <- list("sdev" = sv$d, "rotation" =sv$v,center = cm)

  result$scale <- if(!scale.) FALSE else sds;

  if(retX)
    result$x <- pcscores

  class(result) <- c("pPCA","prcomp")
  return (result)
}



