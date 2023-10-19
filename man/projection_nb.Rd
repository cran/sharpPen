\name{projection_nb}
\alias{projection_nb}
\title{
Projection operator for norm balls.
}
\description{
Compute the projection operator for norm balls. For example,
we construct \deqn{\lambda P_{B_{\| \cdot \|_*}[0,r]}(x/\lambda) = projection_nb(\lambda,r,\| \cdot \|_*,x)}, 
where \eqn{\| \cdot \|_*} can be \eqn{l_{1}}-norm, \eqn{l_{2}}-norm, and 
\eqn{l_{\infty}}-norm.
}
\usage{
projection_nb(
lambda,radius,family=c("norm2","norm1","norminf"),
input)
}

\arguments{
  \item{lambda}{parameter \eqn{\lambda} in the above equation}
  \item{radius}{parameter \eqn{r} in the above equation}
  \item{family}{select the norm ball type, can be \eqn{l_{1}}-norm, \eqn{l_{2}}-norm, and
\eqn{l_{\infty}}-norm.}
  \item{input}{input x in the above equation}
}
\details{
  Take \eqn{x} as input, \eqn{\lambda} and \eqn{r} as parameters.
  Calculate \eqn{\lambda P_{B_{\| \cdot \|_*}[0,r]}(x/\lambda)} 
  for a given norm ball type.
}
\value{
  \item{projection}{\eqn{\lambda P_{B_{\| \cdot \|_*}[0,r]}(x/\lambda)}}
}
\author{
D.Wang and W.J.Braun
}

\examples{
set.seed(1234567)
family <- "norm1"
temp_p1<-rep(10,100)
projection_nb(3,1,family=family,temp_p1)
}
\keyword{projection}
