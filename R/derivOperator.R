derivOperator<-function(penalty,gamma,h, xx,zz,p){
  B<-matrix(0,ncol=length(zz),nrow=length(xx))
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
  getbb<-function(penalty,gamma,h, xx,x,p){
    
    aaFunt<-function(x) {getaa(h,xx,x,p)}
    bb<-rep(0,length(xx))
    if(penalty=="drv1"){
      bb<-numericalDerivative(x,aaFunt,k=1)
    }
    if(penalty=="drv2"){
      bb<-numericalDerivative(x,aaFunt,k=2)
    }
    if(penalty=="drv3"){
      bb<-numericalDerivative(x,aaFunt,k=3)
    }
    if(penalty=="drv4"){
      bb<-numericalDerivative(x,aaFunt,k=4)
    }
    if(penalty=="Exponential"){
      bb<-numericalDerivative(x,aaFunt,k=2)+gamma*numericalDerivative(x,aaFunt,k=1)
    }
    if(penalty=="Periodicity"){
      bb<-numericalDerivative(x,aaFunt,k=4)+gamma*numericalDerivative(x,aaFunt,k=2)
    }
    if(penalty=="Periodicity2"){
    bb<-numericalDerivative(x,aaFunt,k=2)+gamma*numericalDerivative(x,aaFunt,k=0)
    }      
    bb
  }
  
  for(j in 1:length(zz))
    B[,j]<-getbb(penalty,gamma,h,xx,zz[j],p)
  B
}
