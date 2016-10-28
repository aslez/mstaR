print.mstaR <- function(object) {
  cat("\nCall:\n")
  print(object$call)
  cat("\nIndependent variables:\n")
  print(object$coefficients)
  if (!is.null(object$phi)) {
  cat("\nTemporal lag effects:\n")
  print(object$phi)  
  }
  cat("\nSpatial lag effects:\n")
  print(object$rho)
}

print.summary.mstaR <- function(object) {
  cat("\nCall:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients)
  cat("\nLog Likelihood:", formatC(object$log.lik, digits = 5))
  cat("\nAIC:", formatC(object$aic, digits = 5))
  cat("\nSigma:", formatC(object$sigma, digits = 5))
}
