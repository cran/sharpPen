\name{numericalDerivative}
\alias{numericalDerivative}
\title{Numerical Derivative of Smooth Function}
\description{
  Cubic spline interpolation of columns of a matrix for
purpose of computing numerical derivatives at a corresponding
sequence of gridpoints.
}
\usage{
numericalDerivative(x, g, k, delta=.001) 
}
\arguments{
\item{x}{numeric vector}
\item{g}{numeric-valued function of x}
\item{k}{number of derivatives to be computed}
\item{delta}{denominator of Newton quotient approximation}
}
\value{
numeric vector of kth derivative of g(x)
}
\author{W.J. Braun}
\keyword{models}
