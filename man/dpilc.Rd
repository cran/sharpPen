\name{dpilc}
\alias{dpilc}
\title{
Select a Bandwidth for Local Quadratic and Cubic Regression
}
\description{
Use direct plug-in methodology to select the bandwidth
of a local quadratic and local cubic Gaussian kernel regression estimate, 
as an extension of Wand's \code{dpill} function. 
}
\usage{
dpilc(xx, yy, blockmax = 5, divisor = 20, trim = 0.01,
                 proptrun = 0.05, gridsize = 401L, range.x = range(x))
}
\arguments{
\item{xx}{
numeric vector of x data.
Missing values are not accepted.
}
\item{yy}{
numeric vector of y data.
This must be same length as \code{x}, and
missing values are not accepted.
}
\item{blockmax}{
the maximum number of blocks of the data for construction
of an initial parametric estimate. 
}
\item{divisor}{
the value that the sample size is divided by to determine
a lower limit on the number of blocks of the data for
construction of an initial parametric estimate.
}
\item{trim}{
the proportion of the sample trimmed from each end in the
\code{x} direction before application of the plug-in methodology.
}
\item{proptrun}{
the proportion of the range of \code{x} at each end truncated in the
functional estimates.
}
\item{gridsize}{
number of equally-spaced grid points over which the
function is to be estimated.
}
\item{range.x}{
vector containing the minimum and maximum values of \code{x} at which to
compute the estimate.
For density estimation the default is the minimum and maximum data values
with 5\% of the range added to each end.
For regression estimation the default is the minimum and maximum data values.
}

}
\value{
the selected bandwidth.
}
\details{
This function is a local cubic (also quadratic) extension of
the \code{dpill} function of Wand's KernSmooth package.  
The kernel is the standard normal density.
Least squares octic fits over blocks of data are used to 
obtain an initial estimate. As in Wand's implementation
of the Ruppert, Sheather and Wand selector, 
Mallow's \eqn{C_p}{Cp} is used to select
the number of blocks.  An option is available to 
make use of a periodic penalty (with possible trend) 
relating the 4th derivative of the regression function
to a constant (gamma) times the 2nd derivative.  This
avoids the need to calculate the octic fits and reverts
back to the original quartic fits of \code{dpill} with
appropriate adjustments to the estimated functionals
needed in the direct-plug-in bandwidth calculation. This 
code is similar to \code{dpilq} but uses a 6th degree
polyomial approximation instead of an 8th degree polynomial
approximation.
}
\section{Warning}{
If there are severe irregularities (i.e. outliers, sparse regions)
in the \code{x} values then the local polynomial smooths required for the
bandwidth selection algorithm may become degenerate and the function
will crash. Outliers in the \code{y} direction may lead to deterioration
of the quality of the selected bandwidth.
}
\references{
Ruppert, D., Sheather, S. J. and Wand, M. P. (1995).
An effective bandwidth selector for local least squares
regression.
\emph{Journal of the American Statistical Association},
\bold{90}, 1257--1270.

Wand, M. P. and Jones, M. C. (1995).
\emph{Kernel Smoothing.}
Chapman and Hall, London.
}
\seealso{
\code{\link[stats]{ksmooth}}, \code{\link[KernSmooth]{locpoly}}.
}
\examples{
x <- faithful$eruptions
y <- faithful$waiting
plot(x, y)
h <- dpill(x, y)
fit <- locpoly(x, y, bandwidth = h, degree=1)
lines(fit)
h <- dpilc(x, y)
fit <- locpoly(x, y, bandwidth = h, degree=2)
lines(fit, col=3, lwd=2)
fit <- locpoly(x, y, bandwidth = h, degree=3)
lines(fit, col=2, lwd=2)
}
\author{D.Wang and W.J.Braun}
\keyword{models}
