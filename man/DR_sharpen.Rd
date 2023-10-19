\name{DR_sharpen}
\alias{DR_sharpen}
\title{
Shape-Constrained Local Linear Regression via Douglas-Rachford 
}
\description{
Local linear regression is applied to bivariate
data.  The response is `sharpened' or perturbed in a way to render a
curve estimate that satisfies some specified shape constraints. 
}
\usage{
DR_sharpen(
x, y, xgrid=NULL, M=200,  h=NULL, mode=NULL, 
    ratio_1=0.14,ratio_2=0.14,ratio_3=0.14,ratio_4=0.14,
    constraint_1=NULL, constraint_2=NULL, constraint_3=NULL,
    constraint_4=NULL, norm="l2", augmentation=FALSE, maxit = 10^5)
}

\arguments{
  \item{x}{a vector of explanatory variable observations}
  \item{y}{binary vector of responses}
  \item{xgrid}{gridpoints on x-axis where estimates are taken}
  \item{M}{number of equally-spaced gridpoints (if xgrid not specified)}
  \item{h}{bandwidth}
  \item{mode}{the location of the peak on the x-axis, valid in the unimode case}
  \item{ratio_1}{control the first derivative shape constraint gap aroud the peak, 
  valid in the unimode case}
  \item{ratio_2}{control the second derivative shape constraint gap aroud the peak, 
  valid in the unimode case}
  \item{ratio_3}{control the third derivative shape constraint gap aroud the peak, 
  valid in the unimode case}
  \item{ratio_4}{control the fourth derivative shape constraint gap aroud the peak,
  valid in the unimode case}
  \item{constraint_1}{a vector of the first derivative shape constraint}
  \item{constraint_2}{a vector of the second derivative shape constraint}
  \item{constraint_3}{a vector of the third derivative shape constraint}
  \item{constraint_4}{a vector of the fourth derivative shape constraint}
  \item{norm}{the smallest possible distance type: "l2", "l1" or "linf". Default is "l2"}
  \item{augmentation}{data augmentation: "TRUE" or "FALSE", default is "FALSE"}
  \item{maxit}{maximum iterarion number, default is 10^5}
}
\details{
  Data are perturbed the smallest possible L2 or L1 or Linf distance subject to the 
constraint that the local linear estimate satisfies some specified 
shape constraints.
}
\value{
  \item{ysharp}{sharpened responses}
  \item{iteration}{number of iterations the function has been 
  spend for the convergence}
}
\references{
Wang, D. (2022). Penalized and constrained data sharpening methods for kernel regression (Doctoral dissertation, University of British Columbia).
}
\author{
D.Wang and W.J.Braun
}

\examples{
set.seed(1234567)
gam<-4
g <- function(x) (3*sin(x*(gam*pi))+5*cos(x*(gam*pi))+6*x)*x
n<-100
M<-200
noise <- 1
x<-sort(runif(n,0,1))
y<-g(x)+rnorm(n,sd=noise)
z<- seq(min(x)+1/M, max(x)-1/M, length=M) ############xgrid points
h1<-dpill(x,y)
A<-lprOperator(h=h1,x=x,z=z,p=1)
local_fit<-t(A)%*%y
ss_1<-c(sign(numericalDerivative(z,g,k=1)))
DR_sharpen(x=x, y=y, xgrid=z, h=h1, constraint_1=ss_1, norm="linf",maxit =10^3)
}
\keyword{models}
