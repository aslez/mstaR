asvcov <- function(rho, W, b, s2, X) { 
  zvec <- rep(0, length(b))
  
  #sigma
  n <- nrow(W[[1]])
  sigma.list <- lapply(seq_along(W), function(x) sigma.se(W[[x]], rho[x], s2))
  sigma.line <- c(n / 2, unlist(sigma.list), zvec)
  
  #rho
  rho.list <- lapply(seq_along(W), function(x) {
    rho.se(rho[x], rho, W[[x]], W, b, s2, X)
  })
  rho.line <- do.call('rbind', rho.list)
  
  #beta
  beta.list <- lapply(seq_along(W), function(x){
    beta.se(rho[x], W[[x]], s2, b, X)
  })
  beta.line <- cbind(zvec, do.call('cbind', beta.list), s2 * t(X) %*% X)
  
  #result
  vcov.mat <- s2^2 * solve(rbind(sigma.line, rho.line, beta.line))
  dnames <- c('sigma', paste0('Rho.', seq_along(W)), colnames(X))
  rownames(vcov.mat) <- colnames(vcov.mat) <- dnames
  vcov.mat
}

tr <- function(x) sum(diag(x))

sigma.se <- function(w, p, s2) {
  n <- nrow(w)
  Ai <- solve(diag(n) - p * w) 
  B <- Ai %*% w
  s2 * tr(B)
}

rho.se <- function(p, rho, w, W, b, s2, X) {
  n <- nrow(w)
  Ai <- solve(diag(n) - p * w) 
  B1 <- Ai %*% w
  rsein <- lapply(seq_along(W), function(x) {
    rho.se.inner(rho[x], W[[x]], s2, b, X, B1)
  })
  c(s2 * tr(B1), unlist(rsein), s2 * t(X) %*% B1 %*% X %*% b) 
}

rho.se.inner <- function(p, w, s2, b, X, B1) {
  n <- nrow(w)
  Ai <- solve(diag(n) - p * w)
  B2 <- Ai %*% w
  s2^2 * (tr(t(B1) %*% B2) + tr(B1 %*% B2)) + 
    s2 %*% t(b) %*% t(X) %*% t(B1) %*% B2 %*% X %*% b
}

beta.se <- function(p, w, s2, b, X) {
  n <- nrow(w)
  Ai <- solve(diag(n) - p * w) 
  B <- Ai %*% w
  t(s2 * t(b) %*% t(X) %*% t(B) %*% X)
}