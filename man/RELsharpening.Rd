\name{RELsharpening}
\alias{RELsharpening}
\title{Ridge/Enet/LASSO Sharpening via the mean/local polynomial regression with 
       large bandwidth/linear regression.}
\description{
This is a function to shrink responses towards their
mean/estimations of local polynomial regression with large
bandwidth/estimations of linear regression as a form of data sharpening to remove roughness, and reduce the 
bias (when "combine=TRUE"),
prior to use in local polynomial regression.
}
\usage{
RELsharpening(x,y,alpha,type,bigh,hband,combine)
}
\arguments{
\item{x}{numeric vector of equally spaced x data. Missing values are not accepted.}
\item{y}{vector of y data. Missing values are not accepted.}
\item{alpha}{the elasticnet mixing parameter vector, with alpha in [0,1].}
\item{type}{The type of the base line. In total, we have three types: "mean", "big_h", and "linear".} 
\item{bigh}{the kernel bandwidth smoothing parameter.}
\item{hband}{the kernel bandwidth smoothing parameter, which will be used in the residual sharpening method.}
\item{combine}{Should the smoother combined with residual method or not, default=FALSE.}
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
ys<-RELsharpening(x, y,alpha=c(0.2,0.8),"big_h", dpill(x,y)*4, dpill(x,y),combine=TRUE)
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
