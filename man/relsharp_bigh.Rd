\name{relsharp_bigh}
\alias{relsharp_bigh}
\title{Ridge/Enet/LASSO Sharpening via the local polynomial regression with 
       large bandwidth.}
\description{
This is a function to shrink responses towards their
estimations of local polynomial regression with large
bandwidth as a form of data sharpening to remove roughness,
prior to use in local polynomial regression.
}
\usage{
 relsharp_bigh(x, y, alpha, bigh) 
}
\arguments{
\item{x}{numeric vector of equally spaced x data. Missing values are not accepted.}
\item{y}{vector of y data. Missing values are not accepted.}
\item{alpha}{the elasticnet mixing parameter vector, with alpha in [0,1].}
\item{bigh}{the kernel bandwidth smoothing parameter.}
}
\value{
numeric matrix of sharpened responses, with each column
 corresponding to different values of alpha
}
\details{Note that the predictor values are assumed to be equally
spaced.  
}
\author{D.Wang}
\examples{
x<-seq(0,10,length=100)
g <- function(x) sin(x)
y<-g(x)+rnorm(100)
ys<-relsharp_bigh(x, y,alpha=c(0.2,0.8), dpill(x,y)*4)
y.lp2<-locpoly(x,ys[,1],bandwidth=dpill(x,y),degree=1,gridsize=100)
y.lp8<-locpoly(x,ys[,2],bandwidth=dpill(x,y),degree=1,gridsize=100)
y.lp<-locpoly(x,y,bandwidth=dpill(x,y),degree=1,gridsize=100)
curve(g,x,xlim=c(0,10))
lines(y.lp2,col=2)
lines(y.lp8,col=3)
lines(y.lp,col=5)
norm(as.matrix(g(x) - y.lp2$y),type="2")
norm(as.matrix(g(x) - y.lp8$y),type="2")
norm(as.matrix(g(x) - y.lp$y),type="2")
}
\keyword{models}
