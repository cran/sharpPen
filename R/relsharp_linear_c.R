relsharp_linear_c <- function(x, y, alpha, hband) {
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
        bw <- npregbw(formula=newy_l~x,regtype="ll", bwmethod="cv.aic")
        bw$bw <- hband
        model<-npreg(bws = bw,residuals=TRUE)
        y_l_resi<-y+model$resid
        ysharpMat[,i]<-y_l_resi
      }
ysharpMat
}
