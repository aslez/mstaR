summary.mstaR <- function(object) {
  ans <- object["call"]
  if (!is.null(object$phi)) {
    est <- c(object$coefficients, object$phi, object$rho)
    se <- c(object$coef.se, object$phi.se, object$rho.se)
    dnames <- list(c(names(object$coefficients),
                     names(object$phi),
                     names(object$rho)),
                   c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)'))
  }
  else {
    est <- c(object$coefficients, object$rho)
    se <- c(object$coef.se, object$rho.se) 
    dnames <- list(c(names(object$coefficients),
                     names(object$rho)),
                   c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)'))
  }
  zval <- est/se
  ans$rank <- object$rank
  ans$df <- c(object$rank, object$df.residual)
  ans$log.lik <- object$log.lik
  ans$aic <- (2 * object$rank) - (2*ans$log.lik)
  ans$sigma <- sqrt(object$s2)
  names(ans$aic) <- 'AIC'
  ans$lag.count <- length(object$rho)
  ans$coefficients <- cbind(est, se, zval, 2 * pnorm(-abs(zval)))
  dimnames(ans$coefficients) <- dnames
  class(ans) <- 'summary.mstaR'
  ans
}