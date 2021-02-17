\name{getB}
\alias{getB}
\title{Shape Constraint Matrix Construction}
\description{
  Construct a shape constraint matrix at 
  a corresponding sequence of x data and sequence of gridpoint z.
}
\usage{
  getB(penalty,gamma,h, xx,zz,p) 
}
\arguments{
  \item{penalty}{the type of shape constraint, can be "R1", 
    "Roughness", "Exponential" and "Periodicity".}
  \item{gamma}{the shape constraint parameter}
  \item{h}{the kernel bandwidth smoothing parameter.}
  \item{xx}{numeric vector of x data. Missing values are not accepted.}
  \item{zz}{numeric vector of gridpoint z data. Missing values are not accepted.}
  \item{p}{degree of local polynomial used.}
  
}
\value{
    shape constraint matrix 
}
\author{X.J. Hu}
\keyword{models}

