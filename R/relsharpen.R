relsharpen <- function (x, y, h, alpha = c(0, 0.5, 1), p = 2, M = 51) 
{
    if (missing(h)) h <- dpill(x, y)
    constraintpoints <- seq(min(x) + 1/M, max(x) - 1/M, length = M)
    B <- derivOperator(penalty="drv2",gamma=NULL,h, xx = x, zz = constraintpoints, p)
    yp <- t(B) %*% y
    ysharpMat <- matrix(0, nrow = length(x), ncol = length(alpha))
    for (i in 1:length(alpha)) {
        fit <- cv.glmnet(t(B), yp, alpha = alpha[i], intercept = FALSE)
        coef <- coef(fit, s = fit$lambda.min)
        ysharp <- y - coef[-1]
        ysharpMat[, i] <- ysharp
    }
    ysharpMat
}
