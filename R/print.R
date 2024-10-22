#' Print the Output of Principal Component Analysis (pPCA)
#' @rdname print.pPCA
#' @description Prints the output of the \code{pPCA}
#'
#' @param x An object of class \code{pPCA} that contains the results of a partial principal component analysis.
#' @param digits The number of decimal places to use in printing results such as variance explained and PC scores. Defaults to 3.
#' @param \dots Further arguments passed to \code{print} for additional control over the output.
#'
#' @return None.
#' @export

print.pPCA <- function(x, digits = 3, ...) {

  n <- length(x$sdev)

  cat("\nPrincipal Component Analysis (pPCA) Results\n")
  cat("----------------------------------------------------\n")


  cat("Number of components: ", length(x$sdev), "\n")


  if (x$scale[1]) {
    total_variance <- length(x$center)
  } else {
    total_variance <- sum(x$sds^2)
  }


  var_explained <- x$sdev^2


  proportion_total_var <- var_explained / total_variance
  cum_proportion_total_var <- cumsum(proportion_total_var)


  variance_matrix <- cbind(
    Proportion_of_Variance = round(proportion_total_var, digits),
    Cumulative_Variance = round(cum_proportion_total_var, digits)
  )
  rownames(variance_matrix) <- paste0("PC", seq_along(x$sdev))


  cat("\nVariance Explained by PCs:\n")
  print(variance_matrix)


  if (!is.null(x$x)) {
    cat("\nPrincipal Component Scores (Preview):\n")
    top_rows <- head(x$x, n)
    bottom_rows <- tail(x$x, 2)

    print(round(top_rows, digits), ...)

    cat(".\n")
    cat(".\n")
    cat(".\n")

    colnames(bottom_rows) <-  rep("", ncol(bottom_rows))
    print(round(bottom_rows, digits), ...)
  }


  cat("\nCenter:\n")
  center_preview <- c(head(round(x$center, digits), n), " ...", tail(round(x$center, digits), 2))
  print(center_preview, quote = FALSE)


  if (isFALSE(x$scale)) {
    cat("\nNo scaling was applied\n")
  } else {
    cat("\nScale:\n")
    scale_preview <- c(head(round(x$scale, digits), n), " ...", tail(round(x$scale, digits), 2))
    print(scale_preview, quote = FALSE)
  }

  cat("----------------------------------------------------\n")
  invisible(x)
}
