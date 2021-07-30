mvoutlier.CoDa <- 
function (x, quan = 0.75, alpha = 0.025, col.quantile = c(0, 
    0.05, 0.1, 0.5, 0.9, 0.95, 1), symb.pch = c(3, 3, 16, 1, 
    1), symb.cex = c(1.5, 1, 0.5, 1, 1.5), adaptive = TRUE) 
{


orthbasis <- function (D) 
{
    ilrBase <- NULL
    b <- function(D) {
        m <- matrix(1, nrow = D, ncol = D - 1)
        for (i in 1:(ncol(m) - 1)) {
            m[1:i, i + 1] <- 0
            m[i, i] <- -1
        }
        m[i + 1, i + 1] <- -1
        m
    }
    transform <- function(basis) {
        basis <- as.matrix(basis)
        nc <- ncol(basis)
        D <- nrow(basis)
        isPos <- basis > 0
        isNeg <- basis < 0
        nPos <- matrix(1, D, D) %*% isPos
        nNeg <- matrix(1, D, D) %*% isNeg
        basis <- (isPos * nNeg - isNeg * nPos)
        numb <- sapply(1:nc, function(i) {
            1/sqrt(basis[, i] %*% basis[, i])
        })
        numb <- matrix(numb, ncol = nc, nrow = D, byrow = TRUE)
        basis <- basis * numb
        return(basis)
    }
    basis <- b(D)
    basis <- basis * (-1)
    V <- transform(basis)
    ll <- list(V = V, basisv = basis)
    return(ll)
}

gm <- function (x)
{
    exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
}

pivotCoord <- function (x, pivotvar = 1, base = exp(1))
{
    if (dim(x)[2] < 2)
        stop("data must be of dimension greater equal 2")
    if (any(x < 0))
        stop("negative values not allowed")
    norm <- "sqrt((D-i)/(D-i+1))"
    x.ilr <- matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
    D <- ncol(x)
    w <- which(pivotvar != 1:D)
    x <- x[, c(pivotvar, w)]
    for (i in 1:ncol(x.ilr)) {
       x.ilr[, i] <- eval(parse(text = norm)) * log(apply(as.matrix(x[,
                (i + 1):D]), 1, gm)/(x[, i]), base)
    }
    if (is.data.frame(x))
        x.ilr <- data.frame(x.ilr)
    x.ilr <- data.frame(x.ilr)
    if (all(nchar(colnames(x)) > 1)) {
      for (i in 1:(D - 1)) {
          colnames(x.ilr)[i] <- paste(colnames(x)[i], "_",
          paste(substr(colnames(x)[(i + 1):D], 1, 2),
          collapse = "-"), collapse = "", sep = "")
       }
    }
    return(-x.ilr)
}

############################################################################################

    if (!is.matrix(x) && !is.data.frame(x)) 
        stop("x must be matrix or data.frame")
    if (ncol(x) < 3) 
        stop("x must have at least 3 compositional parts")
    Z <- pivotCoord(x)
    #V <- (-orthbasis(ncol(x))$V)
    V <- (orthbasis(ncol(x))$V)
    Zj <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) {
        Zj[, j] <- pivotCoord(cbind(x[, j], x[, -j]))[, 1]
    }
    dimnames(Zj)[[2]] <- names(x)
    rob <- covMcd(Z, alpha = quan)
    if (adaptive) {
        Zarw <- arw(Z, rob$center, rob$cov, alpha = alpha)
        if (Zarw$cn != Inf) {
            alpha1 <- sqrt(c(Zarw$cn, qchisq(c(0.75, 0.5, 0.25), 
                ncol(Z))))
        }
        else {
            alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), 
                ncol(Z)))
        }
        rd2 <- mahalanobis(Z, center = Zarw$m, cov = Zarw$c)
    }
    else {
        cutoff <- qchisq(1 - alpha, ncol(Z))
        rd2 <- mahalanobis(Z, center = rob$center, cov = rob$cov)
        Zarw <- list(m = rob$center, c = rob$cov, cn = cutoff, 
            w = as.numeric(rd2 < cutoff))
        alpha1 <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(Z)))
    }
    rd <- sqrt(rd2)
    covobj <- list(center = Zarw$m, cov = Zarw$c, n.obs = length(rd), 
        mah = rd)
    Z.pca <- suppressWarnings(princomp(Z, covmat = covobj, cor = FALSE))
    pcaclr <- Z.pca
    eval <- eigen(Zarw$c)$values
    pcaclr$sdev <- sqrt(eval)
    pcaclr$scores <- Z.pca$scores
    pcaclr$loadings <- V %*% Z.pca$loadings
    dimnames(pcaclr$loadings)[[1]] <- names(x)
    pcaobj <- list(method = "robust", eigenvalues = eval, princompOutputClr = pcaclr)
    class(pcaobj) <- "pcaCoDa"
    Zcent <- scale(Zj, center = apply(Zj, 2, median), scale = FALSE)
    eucl <- apply(abs(Zcent), 1, median)
    out <- (!Zarw$w)
    lq <- length(col.quantile)
    colcol <- rev(rainbow(lq - 1, start = 0, end = 0.7))[as.integer(cut(eucl, 
        quantile(eucl, col.quantile, labels = 1:(lq - 1))))]
    colbw <- rev(gray(seq(from = 0.1, to = 0.9, length = lq - 
        1)))[as.integer(cut(eucl, quantile(eucl, col.quantile, 
        labels = 1:(lq - 1))))]
    pchvec <- rep(symb.pch[1], nrow(Zj))
    cexvec <- rep(symb.cex[1], nrow(Zj))
    if (length(symb.pch) == 5 & length(symb.cex) == 5) {
        lalpha <- length(alpha1)
        for (j in 1:(lalpha)) {
            pchvec[rd < alpha1[j]] <- symb.pch[j + 1]
            cexvec[rd < alpha1[j]] <- symb.cex[j + 1]
        }
    }
    mvoutlierCoDa <- list(ilrvariables = Zj, outliers = out, 
        pcaobj = pcaobj, colcol = colcol, colbw = colbw, pchvec = pchvec, 
        cexvec = cexvec)
    class(mvoutlierCoDa) <- "mvoutlierCoDa"
    return(mvoutlierCoDa)
}
