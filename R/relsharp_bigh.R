relsharp_bigh <- function(x, y, alpha, bigh) {
    loc_y_h<-locpoly(x, y, degree=1, bandwidth=bigh,gridsize=length(y))
    y_h<-diff(y-loc_y_h$y)
    ysharpMat<-matrix(0,nrow=length(y),ncol=length(alpha))
    Mat_diff<-matrix(0,nrow=length(y)-1,ncol=length(y))
    for (j in 1:length(y)-1){
        Mat_diff[j,j]<--1
        Mat_diff[j,j+1]<-1
    }
    for (i in 1:length(alpha)){
        #ridge, lasso and e-net
        a<-alpha[i]
        fit_h<-cv.glmnet(Mat_diff,y_h,intercept=FALSE,alpha=a)
        gamma_h<-coef(fit_h,fit_h$lambda.min)[-1]
        newy_h<-y-gamma_h
        ysharpMat[,i]<-newy_h
      }
ysharpMat
}
