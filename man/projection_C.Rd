\name{projection_C}
\alias{projection_C}
\title{
Projection operator for rectangle or nonnegative space 
}
\description{
Compute the projection operator for rectangle or nonnegative space. For example,
we construct \deqn{\lambda P_{C}(x/\lambda) = projection_C(\lambda,C,x)}, 
where \eqn{C} can be rectangle or nonnegative space.
}
\usage{
projection_C(
lambda,family=c("rectangle","nonnegative"),
input,bound=c(-1,1))
}

\arguments{
  \item{lambda}{parameter \eqn{\lambda} in the above equation}
  \item{family}{type of \eqn{C}, can be rectangle or nonnegative space}
  \item{input}{input x in the above equation}
  \item{bound}{lower bound and upper bound for rectangle}
}
\details{
  Take \eqn{x} as input, \eqn{\lambda} as parameter.
  Calculate \eqn{\lambda P_{C}(x/\lambda)} for a given \eqn{C} family 
}
\value{
  \item{projection}{\eqn{\lambda P_{C}(x/\lambda)}}
}
\author{
D.Wang and W.J.Braun
}

\examples{
set.seed(1234567)
family <- "nonnegative"
temp_p1<-runif(10,-1,1)
projection_C(0.5,family=family,temp_p1)
}
\keyword{projection}
