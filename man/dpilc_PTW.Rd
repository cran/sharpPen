\name{dpilc_PTW}
\alias{dpilc_PTW}
\title{dpilc_PTW: Local Polynomial Bandwith Estimation with Blockwise Selection and Pointwise Results}
\description{
The function uses raw or trimmed data, applies grid-based binning, and estimates local bandwidth based on the provided parameters.
}
\usage{
dpilc_PTW(xx, yy, blockmax = 5, divisor = 20, trim = 0.01, proptrun = 0.05, 
          gridsize = 401L, range.x = range(x), use_raw_data = FALSE)
}
\arguments{
  \item{xx}{A numeric vector of x-values.}
  \item{yy}{A numeric vector of y-values, corresponding to the x-values in \code{xx}.}
  \item{blockmax}{An integer specifying the maximum number of blocks to be used in the blockwise selection. Default is 5.}
  \item{divisor}{An integer that controls the block size. Default is 20.}
  \item{trim}{A numeric value between 0 and 1 specifying the proportion of data to be trimmed from the extremes. Default is 0.01.}
  \item{proptrun}{A numeric value between 0 and 1 indicating the proportion of data to be excluded from the running procedure. Default is 0.05.}
  \item{gridsize}{An integer specifying the number of grid points to be used. Default is 401L.}
  \item{range.x}{A numeric vector of length 2 indicating the range over which the smoothing is applied. Default is the range of \code{xx}.}
  \item{use_raw_data}{A logical value indicating whether to use the raw data (TRUE) or trimmed data (FALSE) for analysis (default is FALSE).}
}
\details{
  This function provides a point-wise calculation of the functional theta(4,4)  
at each data point \(x_i\). It employs various auxiliary functions for binning the data, calculating local polynomial estimations, and performing necessary cross-validation. The function returns a point-wise estimate of the relevant quantity, which is used for localized analysis of the data distribution.

  The core methodology used in this function is based on the nonparametric regression framework. For detailed information on the theoretical aspects and derivations, refer to Chapter 3 of Fan and Gijbels (1996).
}
\value{
  A list with three elements:
  \item{x}{A numeric vector of x-values.}
  \item{y}{A numeric vector of y-values.}
  \item{h}{A numeric vector of bandwidth values computed for each corresponding x-value.}
}
\references{
  Fan, J., & Gijbels, I. (1996). \emph{Local Polynomial Modelling and its Applications}. Chapman and Hall/CRC.
}
\seealso{
  \code{\link[KernSmooth]{locpoly}}
}
\examples{
  # Example usage:
  x <- rnorm(100)
  y <- rnorm(100)

  # Run the pointwise estimation
  result <- dpilc_PTW(x, y, blockmax = 5, divisor = 20, trim = 0.01, 
                      proptrun = 0.05, gridsize = 40, range.x = range(x), 
                      use_raw_data = TRUE)
  
  # Inspect the result
  plot(result$x, result$h, type = "l", col = "blue", xlab = "X", ylab = "Bandwidth Estimate")
}
\author{
  D.Wang and W.J.Braun
}
\keyword{models}
