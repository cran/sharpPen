\name{data_sharpening}
\alias{data_sharpening}
\title{
Penalized data sharpening for Local Linear, Quadratic and Cubic Regression 
}
\description{
Penalized data sharpening has been proposed as a way to enforce certain constraints on a local polynomial regression estimator. In addition to a bandwidth, a coefficient of the penalty term is also required. We propose propose systematic approaches for choosing these tuning parameters, in part, by considering the optimization problem from the perspective of ridge regression.  
}
\usage{
data_sharpening(xx,yy,zz,p,h=NULL,gammaest=NULL,penalty,lambda=NULL)
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
\item{zz}{
numeric vector of gridpoint z data. Missing values are not accepted.
}
\item{p}{
degree of local polynomial used.
} 
\item{h}{
the kernel bandwidth smoothing parameter. If NULL, this value
will be estimated by function dpill for Local Linear Regression, and
will be estimated by function dpilc for Local Quadratic and Cubic 
Regression.
}
\item{gammaest}{
the shape constraint parameter. Cannot be NULL for Periodic shape constraint. Can be NULL for Exponential shape constraint.
}
\item{penalty}{
the type of shape constraint, can be "Exponential" and "Periodicity".
}
\item{lambda}{
a coefficient of the penalty term, default is NULL.
}  
}
\value{
the sharpened response variable.
}

\examples{
  set.seed(1234567)
  gam<-4
  gamest<-gam
  g <- function(x) 3*sin(x*(gam*pi))+5*cos(x*(gam*pi))+6*x
  sigma<-3
  xx<-seq(0,1,length=100)
  yy<-g(xx)+rnorm(100,sd=sigma)
  zz<-xx
  h1<-dpilc(xx,yy)
  local_fit<-t(lprOperator(h=h1,xx=xx,zz=zz,p=2))\%*\%yy
  y_sharp<-data_sharpening(xx=xx,yy=yy,zz=zz,p=2,gammaest=gamest,penalty="Periodicity")
  sharp_fit<-t(lprOperator(h=h1,xx=xx,zz=zz,p=2))\%*\%y_sharp
  plot(c(min(xx),max(xx)),c(min(yy)-0.5,max(yy)+0.5),type="n",,xlab="x",ylab="y")
  legend("bottomright",legend=c("curve_local","curve_sharpen"), col=c(1,3),bty="n",pch=c("-","-"))
  lines(xx,local_fit)
  lines(xx,sharp_fit,col=3, lwd=2)
  points(xx,yy,col= rgb(0.8,0.2,0.2,0.2))
}
\author{D.Wang and W.J.Braun}
\keyword{nonparametric}
