netBS <- function(mod, R) {
  slag <- mod$slag
  muvec <- c(mod$s2, mod$rho, mod$coefficients)
  if (!is.null(mod$phi)) muvec <- c(muvec, mod$phi)
  samples <- mvrnorm(R, mu = muvec, Sigma = mod$resvar, empirical = TRUE)
  snames <- colnames(samples)
  s2i <- match('s2', snames)
  inti <- match('(Intercept)', snames)
  rhoi <- grep('Rho.', snames)
  bi <- sapply(colnames(mod$X)[-1], function(x) match(x, snames))
  phii <- NULL
  if (!is.null(mod$phi)) phii <- grep('Phi', snames)
  sres <- apply(samples, 1, processXsamples,
                slag = slag, rhoi = rhoi, bi = bi, phii = phii)
  Ai.sd <- apply(simplify2array(lapply(sres, function(x) x$Ai)), 1:2, sd)
  Ai.p <- pnorm(-abs(mod$Ai / Ai.sd))
  W.sd <- apply(simplify2array(lapply(sres, function(x) x$W)), 1:2, sd)
  W.p <- pnorm(-abs(mod$W / W.sd))
  if (!is.null(mod$phi)) {
    S.sd <- apply(simplify2array(lapply(sres, function(x) x$S)), 1:2, sd)
    S.p <- pnorm(-abs(mod$S / S.sd))
  }
  else S.sd <- S.p <- NULL
  list(Ai.sd = Ai.sd, W.sd = W.sd, S.sd = S.sd,
       Ai.p = Ai.p, W.p = W.p, S.p = S.p) 
}

processXsamples <- function(x, slag, rhoi, bi, phii) {
  rho <- x[rhoi]
  W <- Reduce('+', lapply(seq_along(slag), function(x) rho[x] * slag[[x]]))
  n <- nrow(W)
  idt <- diag(n)
  Ai  <- solve(idt - W)
  if (!is.null(phii)) {
    phi <- x[phii]
    S <- solve(idt - W - phi * idt)
  }
  else S <- NULL
  list(W = W, Ai = Ai, S = S)
}
