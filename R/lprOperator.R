lprOperator<-function(h,xx,zz,p){
  A<-matrix(0,ncol=length(zz),nrow=length(xx))
  getX<-function(p,xx,x){	
    X<-matrix(0,ncol=p+1,nrow=length(xx))
    for(k in 1:(p+1)) X[,k]<-(xx-x)^(k-1)
    X
  }
  getWi<-function(h,xx,x){ 
    xxtmp<-(xx-x)/h
    Wi<-ifelse(abs(xxtmp)<0.0001,1,0)
      Wi <-dnorm(xxtmp, sd=1)
    Wi
  }
  getaa<-function(h,xx,x,p){
    X<-getX(p,xx,x)
    W<-diag(getWi(h,xx,x))
    aa<-(solve(t(X)%*%W%*%X)%*%t(X)%*%W)[1,]
    aa
  }
  for(j in 1:length(zz))
    A[,j]<-getaa(h,xx,zz[j],p)
  A
}
