dpilc_PTW <- function (xx, yy, blockmax = 5, divisor = 20, trim = 0.01, proptrun = 0.05, 
                            gridsize = 401L, range.x = range(x),use_raw_data = FALSE) 
{
  blkest46 <- function(x, y, Nval, q) {
    n <- length(x)
    datmat <- cbind(x, y)
    datmat <- datmat[sort.list(datmat[, 1L]), ]
    x <- datmat[, 1L]
    y <- datmat[, 2L]
    qq <- q + 1L
    xj <- rep(0, n)
    yj <- rep(0, n)
    coef <- rep(0, qq)
    Xmat <- matrix(0, n, qq)
    wk <- rep(0, n)
    qraux <- rep(0, qq)
    sigsqe <- 0
    th44e <- 0
    th46e <- 0
    out <- .Fortran("blkest46", as.double(x), as.double(y), 
                    as.integer(n), as.integer(q), as.integer(qq), as.integer(Nval), 
                    as.double(xj), as.double(yj), as.double(coef), as.double(Xmat), 
                    as.double(wk), as.double(qraux), as.double(sigsqe), 
                    as.double(th44e), as.double(th46e), PACKAGE = "sharpPen")
    list(sigsqe = out[[13]], th44e = out[[14]], th46e = out[[15]])
  }
  rlbin <- function(X, Y, gpoints) {
    rlbinsd <- function(X, Y, n, a, b, M, xcnts, ycnts) {
      xcnts <- rep(0, M)
      ycnts <- rep(0, M)
      delta <- (b - a)/(M - 1)
      for (i in 1:n) {
        lxi <- ((X[i] - a)/delta) + 1
        li <- floor(lxi)
        rem <- lxi - li
        if (li >= 1 & li < M) {
          xcnts[li] <- xcnts[li] + (1 - rem)
          xcnts[li + 1] <- xcnts[li + 1] + rem
          ycnts[li] <- ycnts[li] + (1 - rem) * Y[i]
          ycnts[li + 1] <- ycnts[li + 1] + rem * Y[i]
        }
      }
      return(list(xs = xcnts, ys = ycnts))
    }
    n <- length(X)
    M <- length(gpoints)
    a <- gpoints[1L]
    b <- gpoints[M]
    out <- rlbinsd(as.double(X), as.double(Y), as.integer(n), 
                   as.double(a), as.double(b), as.integer(M), double(M), 
                   double(M))
    list(xcounts = out$xs, ycounts = out$ys)
  }
  linbin <- function(X, gpoints) {
    linbinsd <- function(X, n, a, b, M, gcnts) {
      gcnts <- rep(0, M)
      delta <- (b - a)/(M - 1)
      for (i in 1:n) {
        lxi <- ((X[i] - a)/delta) + 1
        li <- floor(lxi)
        rem <- lxi - li
        if (li >= 1 & li < M) {
          gcnts[li] <- gcnts[li] + (1 - rem)
          gcnts[li + 1] <- gcnts[li + 1] + rem
        }
      }
      return(gs = gcnts)
    }
    n <- length(X)
    M <- length(gpoints)
    a <- gpoints[1L]
    b <- gpoints[M]
    out <- linbinsd(as.double(X), as.integer(n), as.double(a), 
                    as.double(b), as.integer(M), double(M))
    return(gcnts = out)
  }
  sdiagsd <- function(xcnts, delta, hdisc, Lvec, indic, midpts, 
                      M, fkap, ipp, ippp, ss, Smat, work, det, ipvt, Sdg) {
    mid <- Lvec[1] + 1
    midpts[1] <- mid
    fkap[mid] <- 1
    for (j in 1:Lvec[1]) {
      fkap[mid + j] <- exp(-(delta * j/hdisc[1])^2/2)
      fkap[mid - j] <- fkap[mid + j]
    }
    for (k in 1:M) {
      if (xcnts[k] != 0) {
        for (j in max(1, k - Lvec[1]):min(M, k + Lvec[1])) {
          if (indic[j] == 1) {
            fac = 1
            ss[j, 1] <- ss[j, 1] + xcnts[k] * fkap[k - 
                                                     j + midpts[1]]
            for (ii in 2:ippp) {
              fac <- fac * delta * (k - j)
              ss[j, ii] <- ss[j, ii] + xcnts[k] * fkap[k - 
                                                         j + midpts[1]] * fac
            }
          }
        }
      }
    }
    for (k in 1:M) {
      for (i in 1:ipp) {
        j <- 1:ipp
        indss <- i + j - 1
        Smat[i, j] <- ss[k, indss]
      }
      Sdg[k] <- solve(Smat)[1, 1]
    }
    return(Sdg)
  }
  sdiag <- function(x, drv = 0L, degree = 3L, bandwidth, range.x) {
    a <- range.x[1L]
    b <- range.x[2L]
    pp <- degree + 1L
    ppp <- 2L * degree + 1L
    tau <- 4
    xcounts <- x
    M <- length(xcounts)
    gpoints <- seq(a, b, length = M)
    delta <- (b - a)/(M - 1L)
    if (length(bandwidth) == 1L) {
      indic <- rep(1, M)
      Q <- 1L
      Lvec <- rep(floor(tau * bandwidth/delta), Q)
      hdisc <- rep(bandwidth, Q)
    }
    else stop("'bandwidth' must be a scalar")
    dimfkap <- 2L * sum(Lvec) + Q
    fkap <- rep(0, dimfkap)
    midpts <- rep(0, Q)
    ss <- matrix(0, M, ppp)
    Smat <- matrix(0, pp, pp)
    work <- rep(0, pp)
    det <- rep(0, 2L)
    ipvt <- rep(0, pp)
    Sdg <- rep(0, M)
    out <- sdiagsd(xcounts, delta, hdisc, Lvec, indic, midpts, 
                   M, fkap, pp, ppp, ss, Smat, work, det, ipvt, Sdg)
    list(x = gpoints, y = out)
  }
  sstdgsd <- function(xcnts, delta, hdisc, Lvec, indic, midpts, 
                      M, fkap, ipp, ippp, ss, uu, Smat, Umat, work, det, ipvt, 
                      SSTd) {
    mid <- Lvec[1] + 1
    midpts[1] <- mid
    fkap[mid] <- 1
    for (j in 1:Lvec[1]) {
      fkap[mid + j] <- exp(-(delta * j/hdisc[1])^2/2)
      fkap[mid - j] <- fkap[mid + j]
    }
    for (k in 1:M) {
      if (xcnts[k] != 0) {
        for (j in max(1, k - Lvec[1]):min(M, k + Lvec[1])) {
          if (indic[j] == 1) {
            fac = 1
            ss[j, 1] <- ss[j, 1] + xcnts[k] * fkap[k - 
                                                     j + midpts[1]]
            uu[j, 1] <- uu[j, 1] + xcnts[k] * fkap[k - 
                                                     j + midpts[1]]^2
            for (ii in 2:ippp) {
              fac <- fac * delta * (k - j)
              ss[j, ii] <- ss[j, ii] + xcnts[k] * fkap[k - 
                                                         j + midpts[1]] * fac
              uu[j, ii] <- uu[j, ii] + xcnts[k] * fkap[k - 
                                                         j + midpts[1]]^2 * fac
            }
          }
        }
      }
    }
    for (k in 1:M) {
      SSTd[k] <- 0
      for (i in 1:ipp) {
        j <- 1:ipp
        indss <- i + j - 1
        Smat[i, j] <- ss[k, indss]
        Umat[i, j] <- uu[k, indss]
      }
      Smat <- solve(Smat)
      for (i in 1:ipp) {
        for (j in 1:ipp) {
          SSTd[k] <- SSTd[k] + Smat[1, i] * Umat[i, j] * 
            Smat[j, 1]
        }
      }
    }
    SSTd
    return(SSTd)
  }
  sstdiag <- function(x, drv = 0L, degree = 3L, bandwidth, 
                      range.x) {
    a <- range.x[1L]
    b <- range.x[2L]
    pp <- degree + 1L
    ppp <- 2L * degree + 1L
    tau <- 4L
    xcounts <- x
    M <- length(xcounts)
    gpoints <- seq(a, b, length = M)
    delta <- (b - a)/(M - 1L)
    if (length(bandwidth) == 1L) {
      indic <- rep(1, M)
      Q <- 1L
      Lvec <- rep(floor(tau * bandwidth/delta), Q)
      hdisc <- rep(bandwidth, Q)
    }
    else stop("'bandwidth' must be a scalar")
    dimfkap <- 2L * sum(Lvec) + Q
    fkap <- rep(0, dimfkap)
    midpts <- rep(0, Q)
    ss <- matrix(0, M, ppp)
    uu <- matrix(0, M, ppp)
    Smat <- matrix(0, pp, pp)
    Umat <- matrix(0, pp, pp)
    work <- rep(0, pp)
    det <- rep(0, 2L)
    ipvt <- rep(0, pp)
    SSTd <- rep(0, M)
    SSTd <- sstdgsd(xcounts, delta, hdisc, Lvec, indic, midpts, 
                    M, fkap, pp, ppp, ss, uu, Smat, Umat, work, det, 
                    ipvt, SSTd)
    list(x = gpoints, y = SSTd)
  }
  cpsd <- function(X, Y, qq, Nmax) {
    RSS <- rep(0, Nmax)
    n <- length(X)
    for (Nval in 1:Nmax) {
      idiv <- floor(n/Nval)
      for (j in 1:Nval) {
        ilow <- (j - 1) * idiv + 1
        iupp <- j * idiv
        if (j == Nval) {
          iupp <- n
        }
        nj <- iupp - ilow + 1
        Xj <- numeric(nj)
        Yj <- numeric(nj)
        Xj[1:nj] <- X[ilow + (1:nj) - 1]
        Yj[1:nj] <- Y[ilow + (1:nj) - 1]
        Xmat <- cbind(rep(1, nj), matrix(0, nrow = nj, 
                                         ncol = qq - 1))
        for (k in 2:qq) {
          Xmat[, k] <- Xj^(k - 1)
        }
        X.QR <- qr(Xmat, complete = TRUE)
        Q2 <- qr.Q(X.QR, complete = TRUE)[, -c(1:qq)]
        a <- t(Q2) %*% Yj
        RSS[Nval] <- sum(a^2) + RSS[Nval]
      }
    }
    Cpvals <- ((n - qq * Nmax) * RSS/RSS[Nmax]) + 2 * qq * 
      (1:Nmax) - n
    return(Cpvals)
  }
  cpblock <- function(X, Y, Nmax, q) {
    n <- length(X)
    datmat <- cbind(X, Y)
    datmat <- datmat[sort.list(datmat[, 1L]), ]
    X <- datmat[, 1L]
    Y <- datmat[, 2L]
    qq <- q + 1L
    RSS <- rep(0, Nmax)
    Xj <- rep(0, n)
    Yj <- rep(0, n)
    coef <- rep(0, qq)
    Xmat <- matrix(0, n, qq)
    Cpvals <- rep(0, Nmax)
    wk <- rep(0, n)
    qraux <- rep(0, qq)
    Cpvec <- cpsd(X, Y, qq, Nmax)
    which.min(Cpvec)
  }
  xy <- cbind(xx, yy)
  xy <- xy[sort.list(xy[, 1L]), ]
  x <- xy[, 1L]
  y <- xy[, 2L]
  if (use_raw_data){
    #    indlow <- floor(trim * length(x)) + 1
    #    indupp <- length(x) - floor(trim * length(x))
    #    x <- x[indlow:indupp]
    #    y <- y[indlow:indupp]
    n <- length(x)
    M <- gridsize
    a <- range.x[1L]
    b <- range.x[2L]
    gpoints <- seq(a, b, length = M)
    out <- rlbin(x, y, gpoints)
    xcounts <- out$xcounts
    ycounts <- out$ycounts
    Nmax <- max(min(floor(n/divisor), blockmax), 1)
    Nval <- cpblock(x, y, Nmax, 6)
    out <- blkest46(x, y, Nval, 6)
    sigsqQ <- out$sigsqe
    th46Q <- out$th46e
    gamseh <- (sigsqQ * (b - a)/(abs(th46Q) * n))
    Rk45 <- 105/(32 * sqrt(pi))
    uk45 <- 360
    if (th46Q < 0) 
      gamseh <- (((44 - 28) * 22.5) * (Rk45/uk45) * gamseh)^(1/11)
    if (th46Q > 0) 
      gamseh <- (((44 + 28) * 22.5) * (Rk45/uk45) * gamseh)^(1/11)
    mddest <- locpoly(x, y, drv = 4L, bandwidth = gamseh, gridsize = length(x),
                      range.x = range.x)$y
    llow <- floor(proptrun * M) + 1
    lupp <- M - floor(proptrun * M)
    th44kn <- sum((mddest[llow:lupp]^2) * xcounts[llow:lupp])/n
    #  th44kn_point <- (mddest[llow:lupp]^2) * x[llow:lupp]
    th44kn_point <- (mddest^2) * x
    C3K0 <- 7881/(32 * 16^2 * sqrt(2 * pi)) + 27/(8 * sqrt(pi)) - 
      4192/(9 * 64 * sqrt(6 * pi))
    C3K <- (8^3 * C3K0)^(1/17)
    lamseh <- C3K * (((sigsqQ^2) * (b - a)/((th44kn * n)^2))^(1/17))
    mest <- locpoly(xcounts, ycounts, degree = 3, bandwidth = lamseh, 
                    range.x = range.x, binned = TRUE)$y
    Sdg <- sdiag(xcounts, bandwidth = lamseh, range.x = range.x)$y
    SSTdg <- sstdiag(xcounts, bandwidth = lamseh, range.x = range.x)$y
    sigsqn <- sum(y^2) - 2 * sum(mest * ycounts) + sum((mest^2) * 
                                                         xcounts)
    sigsqd <- n - 2 * sum(Sdg * xcounts) + sum(SSTdg * xcounts)
    sigsqkn <- sigsqn/sigsqd
    #result <- (27 * sigsqkn * (b - a)/(4 * sqrt(pi) * th44kn * n))^(1/9)
    
    ###### changed for point-wise
    result_fun <- function(x) (27 * sigsqkn * (b - a)/(4 * sqrt(pi) * x * 
                                                         n))^(1/9)
    
    non_zero_indices <- which(th44kn_point != 0)
    x_temp <- x
    y_temp <- y
    th44kn_point_filtered <- th44kn_point[non_zero_indices]
    
    # Apply result_fun to the filtered non-zero th44kn_point
    h <- result_fun(th44kn_point_filtered)
    
    # Create the result list with filtered x, y, and h
    result <- list(
      x = x_temp[non_zero_indices],  # x values corresponding to non-zero th44kn_point
      y = y_temp[non_zero_indices],  # y values corresponding to non-zero th44kn_point
      h = h            # h values calculated from result_fun applied to filtered th44kn_point
    )
  } else {
    indlow <- floor(trim * length(x)) + 1
    indupp <- length(x) - floor(trim * length(x))
    x <- x[indlow:indupp]
    y <- y[indlow:indupp]
    n <- length(x)
    M <- gridsize
    a <- range.x[1L]
    b <- range.x[2L]
    gpoints <- seq(a, b, length = M)
    out <- rlbin(x, y, gpoints)
    xcounts <- out$xcounts
    ycounts <- out$ycounts
    Nmax <- max(min(floor(n/divisor), blockmax), 1)
    Nval <- cpblock(x, y, Nmax, 6)
    out <- blkest46(x, y, Nval, 6)
    sigsqQ <- out$sigsqe
    th46Q <- out$th46e
    gamseh <- (sigsqQ * (b - a)/(abs(th46Q) * n))
    Rk45 <- 105/(32 * sqrt(pi))
    uk45 <- 360
    if (th46Q < 0) 
      gamseh <- (((44 - 28) * 22.5) * (Rk45/uk45) * gamseh)^(1/11)
    if (th46Q > 0) 
      gamseh <- (((44 + 28) * 22.5) * (Rk45/uk45) * gamseh)^(1/11)
    mddest <- locpoly(xcounts, ycounts, drv = 4L, bandwidth = gamseh, 
                      range.x = range.x, binned = TRUE)$y
    llow <- floor(proptrun * M) + 1
    lupp <- M - floor(proptrun * M)
    th44kn <- sum((mddest[llow:lupp]^2) * xcounts[llow:lupp])/n
    th44kn_point <- (mddest[llow:lupp]^2) * xcounts[llow:lupp] # changed for point_wise
    C3K0 <- 7881/(32 * 16^2 * sqrt(2 * pi)) + 27/(8 * sqrt(pi)) - 
      4192/(9 * 64 * sqrt(6 * pi))
    C3K <- (8^3 * C3K0)^(1/17)
    lamseh <- C3K * (((sigsqQ^2) * (b - a)/((th44kn * n)^2))^(1/17))
    mest <- locpoly(xcounts, ycounts, degree = 3, bandwidth = lamseh, 
                    range.x = range.x, binned = TRUE)$y
    Sdg <- sdiag(xcounts, bandwidth = lamseh, range.x = range.x)$y
    SSTdg <- sstdiag(xcounts, bandwidth = lamseh, range.x = range.x)$y
    sigsqn <- sum(y^2) - 2 * sum(mest * ycounts) + sum((mest^2) * 
                                                         xcounts)
    sigsqd <- n - 2 * sum(Sdg * xcounts) + sum(SSTdg * xcounts)
    sigsqkn <- sigsqn/sigsqd
    #result <- (27 * sigsqkn * (b - a)/(4 * sqrt(pi) * th44kn * n))^(1/9)
    
    ###### changed for point-wise
    result_fun <- function(x) (27 * sigsqkn * (b - a)/(4 * sqrt(pi) * x * 
                                                         n))^(1/9)
    
    non_zero_indices <- which(th44kn_point != 0)
    x_temp <- xcounts[llow:lupp]
    y_temp <- ycounts[llow:lupp]
    th44kn_point_filtered <- th44kn_point[non_zero_indices]
    
    # Apply result_fun to the filtered non-zero th44kn_point
    h <- result_fun(th44kn_point_filtered)
    
    # Create the result list with filtered x, y, and h
    result <- list(
      x = x_temp[non_zero_indices],  # x values corresponding to non-zero th44kn_point
      y = y_temp[non_zero_indices],  # y values corresponding to non-zero th44kn_point
      h = h            # h values calculated from result_fun applied to filtered th44kn_point
    )
    
  }
  return(result) 
}