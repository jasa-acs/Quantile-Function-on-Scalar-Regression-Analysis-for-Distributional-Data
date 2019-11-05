
####  datagen.R is a code to generate simulation dataset
### Input
### Argument		Type						Description
#### p		      Vector of length P  in (0,1)	      Probability grids
### Output yields list of length 6  and say qfregfit  as the output object.
#### Argument	 Type			Description
#### dataset[[1]]  P by 4 matrix    ## true Quant  for 4 groups
#### dataset[[2]]  4 by 4 matrix    ## true Moments (mu1,mu2,mu3,mu4) for 4 groups
#### dataset[[3]]  P by N matrix    ## obs quantiles  (Qy(p))
#### dataset[[4]]  N by A matrix    ## obs covariates (X)

datagen <- function(p = signif(seq(0.001, 0.999, length = 1024), 4)) {
  library(sn)
  library(mvtnorm)

  AR_process <- function(rho, sigma) {
    H <- abs(outer(p, p, "-"))
    V <- sigma * rho^H
    np <- nrow(V)
    V[cbind(1:np, 1:np)] <- V[cbind(1:np, 1:np)] * sigma
    return(V)
  }

  quantiles1 <- function(x, probs = c(0.005, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 0.995)) {
    quantile(x, probs)
  }
  quantiles2 <- function(x, probs = c(0.025, 0.975)) {
    quantile(x, probs)
  }
  quantiles <- function(x, probs = p) {
    quantile(x, probs)
  }

  NQ <- qnorm(p, 0, 1)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))

  alphas <- c(0, 0, 0, -4)
  deltas <- alphas / sqrt(1 + alphas^2)

  omega.fit <- (5^2 / (1 - 2 * deltas[4]^2 / pi))^0.5
  xi.fit <- 3 - omega.fit * deltas[4] * sqrt(2 / pi)
  xis <- c(1, 3, 1, xi.fit)
  omegas <- c(5, 5, 6.5, omega.fit)


  z <- seq(-20, 50, by = 0.1)
  pdf1 <- dsn(z, xi = 1, omega = 5, alpha = 0)
  pdf2 <- dsn(z, xi = 3, omega = 5, alpha = 0)
  pdf3 <- dsn(z, xi = 1, omega = 6.5, alpha = 0)
  pdf4 <- dsn(z, xi = xi.fit, omega = omega.fit, alpha = -4)
  ## pdf4 <- dsn(z,   xi=1, omega=5, alpha=   -4 )


  g.mean <- xis + omegas * deltas * sqrt(2 / pi)
  g.vars <- (omegas^2 * (1 - 2 * deltas^2 / pi))^0.5
  g.skew <- (4 - pi) * (deltas * sqrt(2 / pi))^{
    3
  } / (2 * (1 - 2 * deltas^2 / pi)^{
    3 / 2
  })
  g.kurt <- 2 * (pi - 3) * (deltas * sqrt(2 / pi))^{
    4
  } / ((1 - 2 * deltas^2 / pi)^{
    2
  }) + 3


  qdf1 <- qsn(p, xi = 1, omega = 5, alpha = 0)
  qdf2 <- qsn(p, xi = 3, omega = 5, alpha = 0)
  qdf3 <- qsn(p, xi = 1, omega = 6.5, alpha = 0)
  qdf4 <- qsn(p, xi = xi.fit, omega = omega.fit, alpha = -4)



  groups <- 4
  n.each <- 30
  n <- n.each * groups
  set.seed(1313)
  Sample1 <- rsn(n.each * length(p), xi = 1, omega = 5, alpha = 0)
  Sample2 <- rsn(n.each * length(p), xi = 3, omega = 5, alpha = 0)
  Sample3 <- rsn(n.each * length(p), xi = 1, omega = 6.5, alpha = 0)
  Sample4 <- rsn(n.each * length(p), xi = xi.fit, omega = omega.fit, alpha = -4)



  pure_Samples <- matrix(c(Sample1, Sample2, Sample3, Sample4), ncol = n, nrow = length(p))
  PQY <- round(apply(pure_Samples, 2, quantiles), 6)
  AR_1 <- AR_process(rho = 0.9, sigma = 1) ### dim(AR_1)
  set.seed(713418)
  E_p <- rmvnorm(n, mean = rep(0, dim(AR_1)[2]), sigma = AR_1)
  Samples <- pure_Samples + t(E_p)

  QY <- round(apply(Samples, 2, quantiles), 6)


  within <- matrix(rep(1, n.each), ncol = 1, nrow = n.each)
  between <- diag(c(0, rep(1, groups - 1)), ncol = groups, nrow = groups)


  x <- kronecker(between, within)
  x[, 1] <- 1

  true.moments <- cbind(g.mean, g.vars, g.skew, g.kurt)

  output <- list(true = cbind(qdf1, qdf2, qdf3, qdf4), true.moments = true.moments, qy = QY, x = x)
  return(output)
}

###################################################################################################


quantiles1 <- function(x, probs = c(0.005, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 0.995)) {
  quantile(x, probs)
}
quantiles2 <- function(x, probs = c(0.025, 0.975)) {
  quantile(x, probs)
}
quantiles <- function(x, probs = p) {
  quantile(x, probs)
}

AR_process <- function(rho, sigma) {
  H <- abs(outer(p, p, "-"))
  V <- sigma * rho^H
  np <- nrow(V)
  V[cbind(1:np, 1:np)] <- V[cbind(1:np, 1:np)] * sigma
  return(V)
}
