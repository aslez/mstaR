splagll <- function(rho, y, X, W) {
  if (any(rho < -1) || any(rho > 1)) return(NA)
  n <- nrow(W[[1]])
  I <- diag(n)
  pW <- Reduce('+', lapply(seq_along(W), function(x) rho[x] * W[[x]]))
  A <- I - pW
  ldet <- log(det(A))
  M <- I - X %*% solve(t(X) %*% X) %*% t(X)
  SSE <- t(y) %*% t(A) %*% M %*% A %*% y
  s2 <- SSE / n
  ll <- ldet - ((n / 2) * log(2*pi)) - (n / 2) * log(s2) - (1 / (2 * s2)) * SSE    
  ll
}
