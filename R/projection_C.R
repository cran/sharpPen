projection_C<-function(lambda,family=c("rectangle","nonnegative"),input,bound=c(-1,1)){
  x<-input/lambda
  out<-numeric()
  if(family=="rectangle"){
    for (j in 1:length(x)){
      if(x[j]<=bound[1]){
        out[j]<-bound[1]
      }else if (x[j]>=bound[2]){
        out[j]<-bound[2]
      }else{
        out[j]<-x[j]
      }
    }
  }
  if(family=="nonnegative"){
    out<-x
    out[which(x<0)]<-0
  }
  return(lambda*out)
}
