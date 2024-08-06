### Censored skew-normal regression with delayed entry
### Moser A, Clough-Gorr K, Zwahlen M. (2015) Modeling absolute differences in life expectancy with a censored skew-normal regression approach.
###    PeerJ 3:e1162 https://doi.org/10.7717/peerj.1162

### Author: Andre Moser, November 7, 2017

### Main function: censn
### formula: model formula specification, i.e age1~1
### data: data
### failure: failure/death variable, 1 equals non-censored information (death), 0 censored information (person alive); if is.null(failure) failure is replaced by 1
### ltrun: delayed entry information, e.g. variable "age"
### dp: direct parametrisation=T/F

### Postestimation commands: summary.censn, predict.censn

censn <- function(formula, data, failure=NULL, ltrun=NULL, dp=F, subset, weights, na.action, offset, opt.method="Nelder-Mead", model = F, ...) {
  
  require(sn, warn.conflicts = F)

  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "failure", "ltrun", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  n <- length(y)
  failure <- as.vector(model.extract(mf, "failure"))
  w <- as.vector(model.weights(mf))
  ltrun <- as.vector(model.extract(mf, "ltrun"))

  if (is.null(failure)) failure <- rep(1,n)
  
  if (is.null(w)) w <- rep(1,n)
  
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y)) 
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", length(offset), NROW(y)), domain = NA)
  }
  
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coef = if (is.matrix(y)) matrix(, 0, 3) else numeric(), residuals = y, fitted.values = 0 * y, weights = w)
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  
  else {
    x <- model.matrix(mt, mf, contrasts)
    
    max.gamma1 <- 0.5 * (4 - pi) * (2/(pi - 2))^1.5 - (.Machine$double.eps)^(1/4)
    
    qr.x <- qr(x)
    s <- sqrt(sum(qr.resid(qr.x, y)^2)/n)
    gamma1 <- sum(qr.resid(qr.x, y)^3)/(n * s^3)

    if (abs(gamma1) > max.gamma1) 
      gamma1 <- sign(gamma1) * 0.9 * max.gamma1

    cp.param <- c(qr.coef(qr.x, y), scale=s, shape=gamma1)
   
    if (is.null(ltrun)) {
      fit <- optim(cp.param, fn=sn.lik.cen, y=y, X=x, w=w, fail=failure, hessian=T, method=opt.method)
    }
    
    else {
    
    fit <- optim(cp.param, fn=sn.lik.cen.trun, y=y, X=x, w=w, fail=failure, trun=ltrun, hessian=T, method=opt.method)
   
      
    }
  }

  z <- fit

  names(z$par)[names(z$par)=="(Intercept)"] <- "location (mu)"
  names(z$par)[names(z$par)=="scale"] <- "scale (alpha)"
  names(z$par)[names(z$par)=="shape"] <- "skewness (gamma)"
 
  class(z) <- c("censn")
  
  if (dp==T) {
    warning("DP parametrization is used", call.=F)
    z$par[c("location (mu)","scale (alpha)","skewness (gamma)")] <- cp2dp(z$par[c("location (mu)","scale (alpha)","skewness (gamma)")], family="SN")
    names(z$par)[names(z$par)=="location (mu)"] <- "location (xi)"
    names(z$par)[names(z$par)=="scale (alpha)"] <- "scale (sigma)"
    names(z$par)[names(z$par)=="skewness (gamma)"] <- "shape (psi)"
  }
  
  z$coef <- z$par
  z$par <- NULL
  z$loglik <- -z$value/2
  z$value <- NULL
  z$call <- call
  z$na.action <- attr(mf, "na.action")
  z$counts <- NULL
  z$convergence <- NULL
  z$message <- NULL
  z$offset <- offset
  z$terms <- mt
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$x <- x
  z$y <- y
  z$weights <- w
  z$fitted.values <- c(z$x%*%z$coef[1:ncol(z$x)])
  z$residuals <- z$y-z$fitted.values
  z$skewparam <- ifelse(dp==T, "DP", "CP")
  if (model) 
     z$model <- mf
  z
}

### Likelihood function

sn.lik.cen <- function(cp, y, X, fail, w) {
  n <- length(y)
  if (missing(X)) 
    X <- matrix(rep(1, n), n, 1)
  k <- ncol(X)
  
  if (missing(fail)) 
    fail <- rep(1, length(y))
  
  if (missing(w)) 
    w <- rep(1, length(y))
  
  if (any(w < 0)) 
    stop("weights must be non-negative")
  
  if (abs(cp[k + 2]) > 0.9952717) return(Inf)
  
  dp <- cp2dp(cp, "SN")
  xi <- as.vector(X %*% as.matrix(dp[1:k]))
  
  if (dp[k + 1] <= 0) return(NA)

  logl_uncen <- fail*log(dsn(y, xi=xi, omega=dp[k+1], alpha=dp[k+2]))
  logl_cen <- (1-fail)*(log(1-psn(y, xi=xi, omega=dp[k+1], alpha=dp[k+2])))

  return(-2*sum(w*logl_uncen+w*logl_cen))
  
}

sn.lik.cen.trun <- function(cp, y, X, fail, trun, w) {
  n <- length(y)
  if (missing(X)) 
    X <- matrix(rep(1, n), n, 1)
  k <- ncol(X)
  
  if (missing(fail)) 
    fail <- rep(1, length(y))
  
  if (missing(w)) 
    w <- rep(1, length(y))
  
  if (any(w < 0)) 
    stop("weights must be non-negative")
  
  if (abs(cp[k + 2]) > 0.9952717) 
    return(Inf)
  
  dp <- cp2dp(cp, "SN")

  xi <- as.vector(X %*% as.matrix(dp[1:k]))
  if (dp[k + 1] <= 0) 
    return(NA)
  
  logl_uncen <- fail*(log(dsn(y, xi=xi, omega=dp[k+1], alpha=dp[k+2]))-log(1-psn(trun, xi=xi, omega=dp[k+1], alpha=dp[k+2])))
  logl_cen <- (1-fail)*(log(1-psn(y, xi=xi, omega=dp[k+1], alpha=dp[k+2]))-log(1-psn(trun, xi=xi, omega=dp[k+1], alpha=dp[k+2])))

  return(-2*sum(w*logl_uncen+w*logl_cen))
  
}


### Summary after censn
              
summary.censn <- function(object, alpha=0.05) {
  fishermat <- solve(0.5*object$hessian)
  semat <- sqrt(diag(fishermat))
  p.val <- sprintf("%.5f", 2*(1-pnorm(abs(object$coef/semat))))
    
  out <- data.frame(est=object$coef, se=semat, lci=object$coef-qnorm(1-alpha/2)*semat, uci=object$coef+qnorm(1-alpha/2)*semat, z=object$coef/semat, p=p.val)
  
  names(out)[5:6] <- c("z-ratio", "Pr{>|z|}")
  return(out)
}

### Prediction after censn

predict.censn <- function(obj) {
  if (obj$skewparam=="DP") stop("DP parametrisation used")
  return(mod$x%*%mod$coef[1:ncol(mod$x)])
}


