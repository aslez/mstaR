multisar <- function(formula, slag, data, tlag = NULL) {   
  #set up data
  mf <- lm(formula = formula, data = data, method = 'model.frame')
  y <- model.extract(mf, 'response')
  X <- model.matrix(formula, mf)
  if (nrow(X) != nrow(slag[[1]])) {
    stop('The number of rows in the weights matrix does not match the number of rows in the data.')
  }
  n <- nrow(X)
  if (!is.null(tlag)) {
    z <- as.matrix(model.matrix(tlag, data = data)[, -1])
    colnames(z) <- 'Phi'
    X <- cbind(X, z)
  }
  else z <- NULL
  
  #estimate spatial lag parameters
  ll <- maxLik(splagll, 
               start = rep(0, length(slag)),
               y = y, X = X, W = slag)
  rho <- ll$estimate
  rho.loc <- 2:(1 + length(rho))
  converge <- paste0('Code ', ll$code, ': ', ll$message)
  LL <- ll$maximum
  
  #estimate non-spatial parameters
  W <- Reduce('+', lapply(seq_along(slag), function(x) rho[x] * slag[[x]]))
  idt <- diag(nrow(W))
  lm.lag <- lm((y - W %*% y) ~ X - 1)
  b <- coefficients(lm.lag)
  b.loc <- (length(rho) + 2):(1 + length(rho) + length(b))
  rank <- length(rho) + length(b)
  r <- residuals(lm.lag)
  df <- length(r)
  fit <- y - r
  names(r) <- names(fit)
  names(rho) <- paste0('Rho.', 1:length(slag))
  names(b) <- colnames(X)
  names(LL) <- 'Log Likelihood'
  
  #calculate coefficients and asymptotic standard errors
  SSE <- deviance(lm.lag)
  s2 <- SSE / n
  names(s2) <- c('s2')
  vcov.mat <- asvcov(rho, slag, b, s2, X)
  rho.se <- sqrt(diag(vcov.mat))[rho.loc]
  coef.se <- sqrt(diag(vcov.mat))[b.loc]
  
  #update coefficients
  if (!is.null(tlag)) {
    X <- X[, -ncol(X)]
    phi <- b[length(b)]
    phi.se <- coef.se[length(coef.se)]
    b <- b[-length(b)]
    coef.se <- coef.se[-length(coef.se)]   
    S <- solve(idt - W - phi * idt)
  }
  else phi <- phi.se <- S <- NULL
  
  #return
  call<-match.call()
  result <- list(rho = rho, phi = phi, coefficients = b, s2 = s2,  
                 rho.se = rho.se, phi.se = phi.se, coef.se = coef.se, 
                 resvar = vcov.mat, W = W, Ai = solve(idt - W),
                 S = S, X = X, y = y, z = z, slag = slag,
                 log.lik = LL, residuals = r, fitted.values = fit,
                 rank = rank, df.residual = df - rank, 
                 call = call, converge = converge)
  class(result) <- 'mstaR'
  result
}
