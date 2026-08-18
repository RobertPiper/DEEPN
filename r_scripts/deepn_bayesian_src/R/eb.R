# NOT FUNCTIONAL
eb <- function(z, se, lam=0.05) {
  n <- length(z)
  w <- 1/se^2
  m <- weighted.mean(z, w)
  mu <- lam*m + (1-lam)*z
  v <- 1/(w + lam*w)
  mu[is.na(mu)] <- m
  #v[is.na(v)] <- 1/(lam*)
  #list(mu=mu, v=v)
}
