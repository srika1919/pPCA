#' Print the Output of Principal Component Analysis (pPCA)
#' @rdname print.pPCA
#' @description Prints the output of the \code{pPCA}
#'
#' @param x An object of class \code{pPCA} that contains the results of a partial principal component analysis.
#' @param digits The number of decimal places to use in printing results such as variance explained and PC scores. Defaults to 3.
#' @param n Number of rows to display from the start and end of the PC scores. Defaults to 3.
#' @param scale. Logical value indicating whether scaling was applied during the PCA computation. Defaults to \code{TRUE}.
#' @param \dots Further arguments passed to \code{print} for additional control over the output.
#'
#' @return None.
#' @export

print.pPCA <- function(x, digits = 3, ...) {

  n <- length(x$sdev)

  cat("\nPrincipal Component Analysis (pPCA) Results\n")
  cat("----------------------------------------------------\n")

  # Print the number of components
  cat("Number of components: ", length(x$sdev), "\n")

  # Calculate total variance based on scaling
  if (x$scale[1]) {
    total_variance <- length(x$center)  # Total variance for scaled data (number of variables)
  } else {
    total_variance <- sum(x$sds^2)  # Total variance for unscaled data (sum of squared standard deviations)
  }

  # Calculate variance explained by selected PCs
  var_explained <- x$sdev^2  # Variance explained by each selected PC

  # Proportion of variance explained relative to total variance
  proportion_total_var <- var_explained / total_variance
  cum_proportion_total_var <- cumsum(proportion_total_var)  # Cumulative proportion of total variance

  # Create a matrix for variance explained
  variance_matrix <- cbind(
    Proportion_of_Variance = round(proportion_total_var, digits),
    Cumulative_Variance = round(cum_proportion_total_var, digits)
  )
  rownames(variance_matrix) <- paste0("PC", seq_along(x$sdev))

  # Print the variance explained matrix with respect to the total dataset variance
  cat("\nVariance Explained by PCs:\n")
  print(variance_matrix)

  # Print the first n rows and last 2 rows of PC scores with vertical dots
  if (!is.null(x$x)) {
    cat("\nPrincipal Component Scores (Preview):\n")
    top_rows <- head(x$x, n)
    bottom_rows <- tail(x$x, 2)


    # Print the top rows
    print(round(top_rows, digits), ...)

    # Print vertical dots to indicate rows in between
    cat(".\n")
    cat(".\n")
    cat(".\n")

    # Print the bottom rows
    print(round(bottom_rows, digits), ...)
  }

  # Print a preview of the center (column means) without quotes
  cat("\nCenter (column means):\n")
  center_preview <- c(head(round(x$center, digits), n), tail(round(x$center, digits), 2))
  print(center_preview, quote = FALSE)

  # Print a preview of the scale (column standard deviations) without quotes
  if (isFALSE(x$scale)) {
    cat("\nNo scaling was applied (scale. = FALSE)\n")
  } else {
    cat("\nScale (column standard deviations):\n")
    scale_preview <- c(head(round(x$scale, digits), n), tail(round(x$scale, digits), 2))
    print(scale_preview, quote = FALSE)
  }

  cat("----------------------------------------------------\n")
  invisible(x)
}
