relsharp_linear <- function(x, y, alpha) {
    lmy<-lm(y~x)
    linear<-lmy$coefficients[1]+x*lmy$coefficients[2]
    y_l<-diff(y-linear)
    ysharpMat<-matrix(0,nrow=length(y),ncol=length(alpha))
    Mat_diff<-matrix(0,nrow=length(y)-1,ncol=length(y))
    for (j in 1:length(y)-1){
        Mat_diff[j,j]<--1
        Mat_diff[j,j+1]<-1
    }
    for (i in 1:length(alpha)){
        #ridge, lasso and e-net
        a<-alpha[i]
        fit_l<-cv.glmnet(Mat_diff,y_l,intercept=FALSE,alpha=a)
        gamma_l<-coef(fit_l,fit_l$lambda.min)[-1]
        newy_l<-y-gamma_l
        ysharpMat[,i]<-newy_l
      }
ysharpMat
}
