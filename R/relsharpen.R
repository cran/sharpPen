relsharpen <- function (x, y, h, alpha = c(0, 0.5, 1), p = 2, M = 51) 
{
    getB <- function(h, xx, zz, p) {
        B <- matrix(0, ncol = length(zz), nrow = length(xx))
        for (j in 1:length(zz)) {
            B[, j] <- getbb(h, xx, zz[j], p)
        }
        B
    }
    getbb <- function(h, xx, x, p) {
        aaFunt <- function(x) {
            getaa(h, xx, x, p)
        }
        bb <- numericalDerivative(x, aaFunt, k = 2)
        bb
    }
    getaa <- function(h, xx, x, p) {
        X <- getX(p, xx, x)
        W <- diag(getWi(h, xx, x))
        aa <- (ginv(t(X) %*% W %*% X) %*% t(X) %*% W)[1, ]
        aa
    }
    getWi <- function(h, xx, x) {
        xxtmp <- (xx - x)/h
        Wi <- ifelse(abs(xxtmp) < 1e-04, 1, 0)
        Wi <- dnorm(xxtmp, sd = 1)
        Wi
    }
    getX <- function(p, xx, x) {
        X <- matrix(0, ncol = p + 1, nrow = length(xx))
        for (k in 1:(p + 1)) X[, k] <- (xx - x)^(k - 1)
        X
    }
    if (missing(h)) h <- dpill(x, y)
    constraintpoints <- seq(min(x) + 1/M, max(x) - 1/M, length = M)
    B <- getB(h, xx = x, zz = constraintpoints, p)
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
