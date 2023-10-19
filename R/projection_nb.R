projection_nb<-function(lambda,radius,family=c("norm2","norm1","norminf"),input){
  x<-input/lambda
  out<-numeric()
  if(family=="norm2"){
    out<-radius*(x)/max(radius,norm(x,type="f"))
  }
  if (family=="norm1"){
    if(sum(abs(x))<=radius){
      lambda_temp<-0
    }else{
      xs<-sort(abs(x),decreasing=TRUE)
      abssum<-xs[1]-xs[2]
      i<-2
      while(abssum<radius){
        if(i<length(xs)){
          abssum<-abssum+i*(xs[i]-xs[i+1])
          i<-i+1
        }else{
          break
        }
      }
      xss<-xs[1:i]
      lambda_temp<-(sum(xss)-radius)/i
    }
    for (j in 1:length(x)){
      if(abs(x[j])>lambda_temp){
        out[j]<-sign(x[j])*(abs(x[j])-lambda_temp)
      }else{
        out[j]<-0
      }
    }
  }
  if (family=="norminf"){
    out<-x
    if(max(abs(x))<radius){
      out<-x
    }else{
      out[which(x>radius)]<-radius
      out[which(x<(-radius))]<--radius
      
    }
  }
  return(lambda*out)
}
