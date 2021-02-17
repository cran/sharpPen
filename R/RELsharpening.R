RELsharpening<-function(x,y,alpha=c(0,0.5,1),type,bigh,hband,combine=FALSE){
  ysharpMat<-matrix(0,nrow=length(y),ncol=length(alpha))
  Mat_diff<-matrix(0,nrow=length(y)-1,ncol=length(y))
  for (j in 1:length(y)-1){
  Mat_diff[j,j]<--1
  Mat_diff[j,j+1]<-1
  }
  
  if(missing(combine)){
    if (type=="mean") ysharpMat <- relsharp_mean(y, alpha)
    if (type=="big_h") ysharpMat <- relsharp_bigh(x, y, alpha, bigh)
    if (type=="linear") ysharpMat <- relsharp_linear(x, y, alpha)
  }   
  if(combine==TRUE)  {
    if (type=="mean") ysharpMat <- relsharp_mean_c(x, y, alpha, hband)
    if (type=="big_h")ysharpMat <- relsharp_bigh_c(x, y, alpha, bigh, hband)
    if (type=="linear") ysharpMat <- relsharp_mean_c(x, y, alpha, hband)
    
  }
  ysharpMat
}
