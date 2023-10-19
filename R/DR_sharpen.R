DR_sharpen<-function(
    x, y, xgrid=NULL, M=200,  h=NULL, mode=NULL, 
    ratio_1=0.14,ratio_2=0.14,ratio_3=0.14,ratio_4=0.14,
    constraint_1=NULL, constraint_2=NULL, constraint_3=NULL,
    constraint_4=NULL, norm="l2", augmentation=FALSE, maxit = 10^5){

   
bandwidth<-function(x,y){
    if ((is.numeric(try(dpill(x,y), silent = TRUE))) & (!is.na(try(dpill(x,y), silent = TRUE)))){
      h <-dpill(x,y)
    }else{
      z<-seq(min(x), max(x), length=length(x))
      h<-dpill(z,y)
    }
  return(h)
}

  
  if (augmentation){
    gap <- max(diff(x))
    gap_thresh <- thumbBw(x,y,deg=1,kernel=gaussK)
    while (gap > gap_thresh) {
      gap <- max(diff(x))
      index <- which(diff(x)>(gap/2))
      xnew <- apply(cbind(x[index], x[index+1]), 1, mean)
      ynew <- (y[index] + y[index+1])/2
      xx <- c(x, xnew)
      yy <- c(y, ynew)
      y <- yy[order(xx)]
      x <- xx[order(xx)]
    }
  }
  
  if (is.null(h)){
    h <-bandwidth(dpill,x,y)
  }
  
  if (is.null(xgrid)){
    xgrid <- seq(min(x)+1/M, max(x)-1/M, length=M)
  }
  
  M <- length(xgrid)
  n <- length(y)
  
  
  B1 <- derivOperator(penalty="drv1",gamma=NULL,h,xx=x,zz=xgrid,p=1)
  B2 <- derivOperator(penalty="drv2",gamma=NULL,h,xx=x,zz=xgrid,p=1)
  B3 <- derivOperator(penalty="drv3",gamma=NULL,h,xx=x,zz=xgrid,p=1)
  B4 <- derivOperator(penalty="drv4",gamma=NULL,h,xx=x,zz=xgrid,p=1)
  
  alpha <- 2
  
  if (is.null(constraint_1)){
    if (is.null(mode)){
      constraint_1 <- rep(0,M)
    }else{
      peak1 <- max(which(xgrid <= mode))
      constraint_1 <- c(rep(1,peak1-ceiling(n*ratio_1)),rep(0,ceiling(n*ratio_1)*2),rep(-1,M-peak1-ceiling(n*ratio_1)))
    }
  }
  
  
  if (is.null(constraint_2)){
    if (is.null(mode)){
      constraint_2 <- rep(0,M)
    }else{
      peak2 <- max(which(xgrid <= mode))
      constraint_2 <- c(rep(1,peak2-ceiling(n*ratio_2)),rep(0,ceiling(n*ratio_2)*2),rep(-1,M-peak2-ceiling(n*ratio_2)))
    }
  }
  
  if (is.null(constraint_3)){
    if (is.null(mode)){
      constraint_3 <- rep(0,M)   
    }else{
      peak3 <- max(which(xgrid <= mode))
      constraint_3 <- c(rep(1,peak3-ceiling(n*ratio_3)),rep(0,ceiling(n*ratio_3)*2),rep(-1,M-peak3-ceiling(n*ratio_3)))
    }
  }
  
  if (is.null(constraint_4)){
    if (is.null(mode)){
      constraint_4 <- rep(0,M)   
    }else{
      peak4 <- max(which(xgrid <= mode))
      constraint_4 <- c(rep(1,peak4-ceiling(n*ratio_4)),rep(0,ceiling(n*ratio_4)*2),rep(-1,M-peak4-ceiling(n*ratio_4)))
    }
  }
  
  k<-1
  x_dr<-rep(0,n)
  xold<-x_dr-1
  estimation<-rep(0,n)
  
  L1<-diag(constraint_1)%*%t(B1)
  r1<--L1%*%y
  v1<-L1%*%x_dr-r1
  sigma1<-0.5
  
  L2<-diag(constraint_2)%*%t(B2)
  r2<--L2%*%y
  v2<-L2%*%x_dr-r2
  sigma2<-0.5
  
  L3<-diag(constraint_3)%*%t(B3)
  r3<--L3%*%y
  v3<-L3%*%x_dr-r3
  sigma3<-0.5
  
  L4<-diag(constraint_4)%*%t(B4)
  r4<--L4%*%y
  v4<-L4%*%x_dr-r4
  sigma4<-0.5
  
  tau<-0.8/(sigma1*norm(L1,type="f")^2+sigma2*norm(L2,type="f")^2+sigma3*norm(L3,type="f")^2+sigma4*norm(L4,type="f")^2)
  
  if (norm=="l1"){
    family = "norminf"
  }else if (norm=="l2"){
    family = "norm2"
  }else if (norm=="linf"){
    family = "norm1"
  }
  
  
  while(k<=maxit){
    temp_p1<-x_dr-tau*t(L1)%*%v1-tau*t(L2)%*%v2-tau*t(L3)%*%v3-tau*t(L4)%*%v4
    
    p_1n<-temp_p1-projection_nb(tau,1,family=family,temp_p1)
    w1<-alpha*p_1n-x_dr
    xold<-x_dr
    x_dr<-p_1n
    
    temp_p21<-v1+sigma1*L1%*%w1-sigma1*r1
    p21<-temp_p21-projection_C(sigma1,family="nonnegative",temp_p21)
    v1<-p21
    
    temp_p22<-v2+sigma2*L2%*%w1-sigma2*r2
    p22<-temp_p22-projection_C(sigma2,family="nonnegative",temp_p22)
    v2<-p22
    
    temp_p23<-v3+sigma3*L3%*%w1-sigma3*r3
    p23<-temp_p23-projection_C(sigma3,family="nonnegative",temp_p23)
    v3<-p23
    
    temp_p24<-v4+sigma4*L4%*%w1-sigma4*r4
    p24<-temp_p24-projection_C(sigma4,family="nonnegative",temp_p24)
    v4<-p24
    
    k=k+1
    estimation<-x_dr
  }
  ds_result <- estimation+y
  list(ysharp=ds_result,iteration=k)
}

