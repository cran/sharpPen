relsharp_mean_c <- function(x, y, alpha, hband) {
    y_m<-diff(y)
    ysharpMat<-matrix(0,nrow=length(y),ncol=length(alpha))
    Mat_diff<-matrix(0,nrow=length(y)-1,ncol=length(y))
    for (j in 1:length(y)-1){
        Mat_diff[j,j]<--1
        Mat_diff[j,j+1]<-1
    }
for (i in 1:length(alpha)){
        #ridge, lasso and e-net
        a<-alpha[i]
        fit_m<-cv.glmnet(Mat_diff,y_m,intercept=FALSE,alpha=a)
        gamma_m<-coef(fit_m,fit_m$lambda.min)[-1]
        newy_m<-y-gamma_m
        bw <- npregbw(formula=newy_m~x,regtype="ll", bwmethod="cv.aic")
        bw$bw <- hband
        model<-npreg(bws = bw,residuals=TRUE)
        y_m_resi<-y+model$resid
        ysharpMat[,i]<-y_m_resi
      }
ysharpMat
}
