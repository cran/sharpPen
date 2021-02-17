\name{getA}
\alias{getA}
\title{Local Polynomial Estimator Matrix Construction}
\description{
  Construct a matrix based on the local polynomial estimation at 
  a corresponding sequence of x data and sequence of gridpoint z.
}
\usage{
  getA(h,xx,zz,p)
}
\arguments{
  \item{h}{the kernel bandwidth smoothing parameter.}
  \item{xx}{numeric vector of x data. Missing values are not accepted.}
  \item{zz}{numeric vector of gridpoint z data. Missing values are not accepted.}
  \item{p}{degree of local polynomial used.}
}
\value{
  local polynomial estimator matrix 
}
\author{X.J. Hu}
\keyword{models}

