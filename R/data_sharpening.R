data_sharpening<-function(xx,yy,zz,p,h=NULL,gammaest=NULL,penalty,lambda=NULL){
  n<-length(xx)
  if(penalty=="Periodicity"){
    if(is.null(h)){
      if(p==1){
        h<-dpill(xx,yy)
      }else if(p==2 | p==3){
        h<-dpilc(xx,yy)
      }
    }
    if(is.null(gammaest)){
      stop("'gammaest' cannot be null for periodic penalty, please enter 
            the number of periods this data has.")
    }
    if(is.null(lambda)){
      B1<-derivOperator(penalty="Periodicity",gamma=(gammaest*pi)^2,h,xx=xx,zz=zz,p)
      yp1<-t(B1)%*%yy
      sdy1<-sqrt(var(yp1)*(length(yp1)-1)/length(yp1))[1,1]
      lambda_max <- max(abs(colSums(as.matrix(t(B1))*as.vector(yp1))))/(dim(t(B1))[1]*0.001)
      a<-0.5
      rr<-c(0:50)
      lamseq<-c(lambda_max*a^rr)
      X1<-t(apply(t(B1), MARGIN = 2, FUN = function(X) X/sqrt(sum(diag(X%*%t(X))))))
      kappa_seq<-numeric()
      data_matrix<-matrix(NA,nrow=length(yy),ncol=length(lamseq))
      lamseq_inv<-numeric()
      for(ii in 1:length(lamseq)){
        lamb3<-lamseq[ii]
        lam3<-1/(lamb3*length(yp1))
        lamseq_inv[ii]<-lam3
        if(kappa(1/lam3*diag(n)+B1%*%t(B1),exact=TRUE)<=10^14){
          kappa_seq[ii]<-kappa((1/lam3*diag(n)+X1%*%t(X1)),exact=TRUE)
        }
      }
      lambda<-lamseq_inv[max(which(kappa_seq<=1.000001))]
    }
    y_sharp<-solve(diag(n)+lambda*B1%*%t(B1))%*%yy
    }
  if(penalty=="Exponential"){
    if(is.null(h)){
      if(p==1){
        h<-dpill(xx,yy)
      }else if(p==2 | p==3){
        h<-dpilc(xx,yy)
      }
    }
    if(is.null(gammaest)){
      xy<-as.data.frame(cbind(xx,yy))
      colnames(xy)<-c("x","y")
      fm0<-lm(log(y)~x,xy)
      fm <- nls(y ~ alpha*exp(beta*x)+theta, xy, start = list(alpha=coef(fm0)[1], beta=coef(fm0)[2],theta=min(yy)))
      gammaest<-coef(fm)[2]
    }
    if(is.null(lambda)){
      B1<-derivOperator(penalty="Exponential",gamma=-gammaest,h,xx=xx,zz=zz,p)
      yp1<-t(B1)%*%yy
      sdy1<-sqrt(var(yp1)*(length(yp1)-1)/length(yp1))[1,1]
      lambda_max <- max(abs(colSums(as.matrix(t(B1))*as.vector(yp1))))/(dim(t(B1))[1]*0.001)
      a<-0.5
      rr<-c(0:50)
      lamseq<-c(lambda_max*a^rr)
      X1<-t(apply(t(B1), MARGIN = 2, FUN = function(X) X/sqrt(sum(diag(X%*%t(X))))))
      kappa_seq<-numeric()
      data_matrix<-matrix(NA,nrow=length(yy),ncol=length(lamseq))
      lamseq_inv<-numeric()
      for(ii in 1:length(lamseq)){
        lamb3<-lamseq[ii]
        lam3<-1/(lamb3*length(yp1))
        lamseq_inv[ii]<-lam3
        if(kappa(1/lam3*diag(n)+B1%*%t(B1),exact=TRUE)<=10^14){
          kappa_seq[ii]<-kappa((1/lam3*diag(n)+X1%*%t(X1)),exact=TRUE)
        }
      }
      lambda<-lamseq_inv[max(which(kappa_seq<=5))]
    }
    y_sharp<-solve(diag(n)+lambda*B1%*%t(B1))%*%yy
  }
  return(y_sharp)
  }

