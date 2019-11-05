TruncPCA <- function(input.mat, length.train = NULL, max.K = 20) {
  sizes <- dim(input.mat)
  P <- sizes[1]
  N <- sizes[2]
  ## all.p = seq(P)
  all.n <- seq(N)

  set.seed(12)
  if (is.null(length.train) == TRUE) {
    n <- N
    train.set.n <- all.n
    test.set.n <- NULL
  }
  if (is.null(length.train) != TRUE) {
    n <- length.train
    train.set.n <- sort(sample.int(N, n))
    test.set.n <- all.n[-train.set.n]
  }
  ## if( is.null(n.grid) == TRUE){  p=P;  train.set.p = all.p;    test.set.p = NULL;  }
  ## if( is.null(n.grid) != TRUE){ p=n.grid;  train.set.p = sort( sample.int( P, p ) ); test.set.p =all.p[-train.set.p]; }


  X.train <- t(input.mat[, train.set.n ]) ## dim(X.train)
  X.test <- t(input.mat[, test.set.n  ]) ## dim(X.test)
  n.train <- dim(X.train)[1]
  n.test <- dim(X.test)[1]

  Mu <- apply(X.train, 2, mean)
  Mu.train <- matrix(rep(Mu, each = n.train), nrow = n.train, ncol = P)
  Mu.test <- matrix(rep(Mu, each = n.test), nrow = n.test, ncol = P)

  Res.train <- X.train - Mu.train
  Res.test <- X.test - Mu.test


  Cov <- t(Res.train) %*% (Res.train) / (n.train - 1) ## dim(Cov)
  svd.Cov <- svd(Cov)
  D2 <- svd.Cov$d
  V <- svd.Cov$v

  cci.train <- matrix(NA, nrow = n.train, ncol = max.K)
  cci.test <- matrix(NA, nrow = n.test, ncol = max.K)
  mse.train <- matrix(NA, nrow = n.train, ncol = max.K)
  mse.test <- matrix(NA, nrow = n.test, ncol = max.K)

  for (i in 1:max.K) { ##  i=1
    phi <- as.matrix(V[, 1:i])
    Proj_phi <- t(Mu.train) + phi %*% solve(t(phi) %*% phi) %*% t(phi) %*% t(Res.train)
    cci.train[, i] <- agree.cci(t(X.train), Proj_phi)
    mse.train[, i] <- apply(abs(t(X.train) - Proj_phi), 2, function(x) sum(x^2)) / P

    if (n.test == 0) {
      cci.test <- NULL
      mse.test <- NULL
    }
    if (n.test != 0) {
      Proj_phi.test <- t(Mu.test) + phi %*% solve(t(phi) %*% phi) %*% t(phi) %*% t(Res.test)
      cci.test[, i] <- agree.cci(t(X.test), Proj_phi.test)
      mse.test[, i] <- apply(abs(t(X.test) - Proj_phi.test), 2, function(x) sum(x^2)) / P
    }
  }
  outputs <- list(V[, 1:max.K], D2[1:max.K], Mu, cci.train, mse.train, cci.test, mse.test)
  names(outputs) <- list("v", "d", "mu", "cci.fit", "mse.fit", "cci.pred", "mse.pred")
  return(outputs)
}

##############################################################################################################
##   Q_center'Q_center = VNDV'
##   u = [ coef_mu , Q_centerV(ND)^-0.5]
##   B' = [ mu , V(ND)^0.5]
##   uB = mu x alpha_mu  +  Q_center x V V'


##   PCscore  =  Ui1
##   BackTransfor =  cbind( Qy_mean , Phi1%*%solve(scale1))

MCMC_PCR <- function(X, PCscore, BackTransfor, n.iter = 2000, burn = 200, ep1 = 0.0064) {
  Sys.time() -> start
  px <- dim(X)[2]
  N <- dim(X)[1]
  Px <- dim(X)[2]

  sd_l2_lse <- PCscore
  K_lse <- dim(sd_l2_lse)[2]
  sPhi_lse <- BackTransfor

  B00_lse <- solve(t(X) %*% X) %*% t(X) %*% sd_l2_lse


  nu0 <- ep1
  ## s20 = 1500
  ## n.iter = 2000

  MCMC_BETA <- matrix(0, ncol = Px * K_lse, nrow = n.iter)
  MCMC_OMEGA <- matrix(0, ncol = K_lse, nrow = n.iter)

  Y_lse <- sd_l2_lse

  for (j in 1:K_lse) { ## j=1
    set.seed(5000 + j)
    d_j <- as.matrix(Y_lse[, j])
    Hg <- X %*% solve(t(X) %*% X) %*% t(X)
    SSRg <- t(d_j) %*% (diag(1, nrow = N) - Hg) %*% d_j
    s2 <- 1 / rgamma(n.iter, (nu0 + N) / 2, (nu0 + SSRg) / 2)
    MCMC_OMEGA[, j] <- s2
    Vb <- solve(t(X) %*% X)
    Eb <- Vb %*% t(X) %*% d_j
    E <- matrix(rnorm(n.iter * Px, 0, sqrt(s2)), n.iter, Px)
    beta <- t(t(E %*% chol(Vb)) + c(Eb))
    s.posit <- (j - 1) * Px + 1
    e.posit <- (j - 1) * Px + Px
    MCMC_BETA[, (s.posit:e.posit) ] <- beta
  }

  POST_BETA_UI <- MCMC_BETA[-c(1:burn), ]
  POST_OMEGA_UI <- MCMC_OMEGA[-c(1:burn), ]
  # QUAN_BETA_UI =  apply(   POST_BETA_UI ,  2 ,  quantiles1)
  # QUAN_OMEGA_UI = apply(  POST_OMEGA_UI ,  2 ,  quantiles1)
  # MEAN_BETA_UI =   apply( POST_BETA_UI, 2 , mean )
  # MEAN_OMEGA_UI =  apply( POST_OMEGA_UI, 2 , mean )
  Stage_Est_1 <- (Sys.time() - start)
  outputs <- list(POST_BETA_UI, Stage_Est_1)
  names(outputs) <- list("postbeta", "times")
  return(outputs)
}

##############################################################################################################



########################################
LCCC.PCR.FIT <- function(Y, B = fit.pca$v, remain.basis) {
  Sys.time() -> start
  Y <- as.matrix(Y)
  n <- ncol(Y)
  active.set <- (remain.basis < min(1000, n))
  feasible.long <- sum(active.set)
  feasible.idx <- seq(length(remain.basis))[ active.set ]

  max.long <- nrow(B)
  Values <- matrix(NA, nrow = n, ncol = feasible.long)

  pcsr_be_0 <- apply(Qy, 1, mean)
  pcsr_be_0_emp <- t(rep(1, n))
  pcsr_fit_0 <- pcsr_be_0 %*% pcsr_be_0_emp

  for (j in 1:feasible.long) { ###  j =2

    pcsr_be <- cbind(apply(Qy, 1, mean), B[, 1:remain.basis[ active.set ][j] ])
    pcsr_be_emp <- rbind(t(rep(1, n)), t(B[, 1:remain.basis[ active.set ][j] ]) %*% (Y - pcsr_fit_0))
    pcsr_fit_v <- pcsr_be %*% pcsr_be_emp
    Values[, j] <- agree.cci(Y, pcsr_fit_v)
  }
  Stage_Est_1 <- (Sys.time() - start)
  outputs <- list(Values, remain.basis[ active.set ], Stage_Est_1)
  return(outputs)
}
#######################################################################

## LCCC.WAVE.FIT --  computes the concordance correlation index for the given basis after denoising process
##  input
## 	global.list: leave-one-out list (it may be output of IncidenceVec function)
## 	B: pre-basis function (Beta CDFs)
## 	remain.counts:  vector containing c value in paper  (it may be output[[2]] of CountBasis)
## 	remain.basis :  vector containing k(c) value in paper  (it may be output[[3]] of CountBasis)
## 	Y.list:  raw.dataset

##  output
## 	list[[1]]: predicted values ( length(probability grids) times n times length(c)  array output ) based on leaveout.list



LCCC.WAVE.FIT <- function(global.list, B = BETA_BASE_TOTAL_2, remain.counts, remain.basis, Y.list) {
  Sys.time() -> start
  active.set <- (remain.basis < 1000)
  feasible.long <- sum(active.set)
  feasible.idx <- seq(length(global.list))[ active.set ]

  max.long <- max(unlist(lapply(Y.list, length)))
  Values <- array(NA, c(max.long, n, feasible.long))


  for (j in 1:feasible.long) { ###  j =21

    remain.counts[active.set ][j]

    REDUCED_BASE_TEMP <- BETA_BASE_TOTAL_2[, global.list[[  feasible.idx[j]  ]][-1]  ] ## dim(REDUCED_BASE_TEMP)

    if (is.null(ncol(REDUCED_BASE_TEMP)) == TRUE) {
      REDUCED_BASE_TEMP <- as.matrix(REDUCED_BASE_TEMP)
    }
    Gram_TEMP <- gramSchmidt(REDUCED_BASE_TEMP, tol = .Machine$double.eps^0.5) ## dim(ORTHON_BASE_TEMP)
    ORTHON_BASE_TEMP0 <- Gram_TEMP$Q
    set.seed(123)
    ORTHON_BASE_TEMP <- ORTHON_BASE_TEMP0

    if (ncol(ORTHON_BASE_TEMP) > 1) {
      fitwavelet <- cbind(ORTHON_BASE_TEMP[, 1], WaveletDenose(ORTHON_BASE_TEMP[, -1], filter.number = 2, family = "DaubExPhase")[, , 2])
      fitwavelet_re <- centering.function(fitwavelet, scale = TRUE)
    }
    if (ncol(ORTHON_BASE_TEMP) == 1) {
      fitwavelet_re <- as.matrix(ORTHON_BASE_TEMP)
    }


    for (i in 1:n) { ### i = 1
      y <- as.vector(Y.list[[i]])
      y.long <- length(y)
      quan_be_0 <- rep(1, y.long)
      quan_be_0_emp <- solve(t(quan_be_0) %*% quan_be_0) %*% t(quan_be_0) %*% y
      quan_fit_0 <- quan_be_0 %*% quan_be_0_emp

      fitwavelet_be <- cbind(quan_be_0, fitwavelet_re)
      # fitwavelet_be_emp =   solve(t(fitwavelet_be )%*%fitwavelet_be )%*%t(fitwavelet_be )%*%y
      fitwavelet_fit_v <- try(fitwavelet_be %*% ginv(t(fitwavelet_be) %*% fitwavelet_be, tol = sqrt(.Machine$double.eps)) %*% t(fitwavelet_be) %*% y)

      Values[(1:y.long), i, j] <- fitwavelet_fit_v
    }
  }
  Stage_Est_1 <- (Sys.time() - start)
  outputs <- list(Values, Stage_Est_1)
  names(outputs) <- list("Values", "times")
  return(outputs)
}

#######################################################################

## MCMC_QR --  computes tau-th quantile(L1)-regression for each p point.
##  input
## 	X: covariates
## 	Y: empirical coefficients
## 	n.iter: # of MCMC iteration
## 	burn:   # of burning
## 	ep1:    prior for nu in the paper (inverse gamma dist)
## 	tau:    set tau in (0,1) for the tau-th quantile regression
##  output
## 	list[[1]]: MCMC samples
## 	list[[2]]: running time


MCMC_QR <- function(X, y, n.iter, burn, ep1, tau) {
  library("quantreg")

  Sys.time() -> start
  px <- dim(X)[2]
  N <- dim(X)[1]
  Px <- dim(X)[2]

  pgird.long <- dim(y)[1]

  nu0 <- ep1
  ## s20 = 1500
  ## n.iter = 2000

  MCMC_BETA <- matrix(0, ncol = Px * pgird.long, nrow = n.iter)
  MCMC_OMEGA <- matrix(0, ncol = pgird.long, nrow = n.iter)

  for (it in 1:n.iter) { ## it = 1

    for (j in 1:pgird.long) { ## j=1
      set.seed(n.iter * (j - 1) + it)
      d_j <- as.matrix(y[j, ])
      rqfit <- rq(d_j ~ -1 + X)
      Hinv <- summary(rqfit, cov = T, se = "ker")$Hinv
      J <- summary(rqfit, cov = T, se = "ker")$J
      cov <- tau * (1 - tau) * Hinv %*% J %*% Hinv
      SSRg <- sum(rqfit$residuals^2)
      s2 <- 1 / rgamma(1, (nu0 + N) / 2, (nu0 + SSRg) / 2)
      MCMC_OMEGA[it, j] <- s2

      Vb <- cov
      Eb <- rqfit$coefficients
      E <- matrix(rnorm(1 * Px, 0, sqrt(s2)), 1, Px)
      beta <- t(t(E %*% chol(Vb)) + c(Eb))
      s.posit <- (j - 1) * Px + 1
      e.posit <- (j - 1) * Px + Px
      MCMC_BETA[ it, (s.posit:e.posit) ] <- beta
    }
  }
  POST_BETA_UI <- MCMC_BETA[-c(1:burn), ]
  POST_OMEGA_UI <- MCMC_OMEGA[-c(1:burn), ]
  # QUAN_BETA_UI =  apply(   POST_BETA_UI ,  2 ,  quantiles1)
  # QUAN_OMEGA_UI = apply(  POST_OMEGA_UI ,  2 ,  quantiles1)
  # MEAN_BETA_UI =   apply( POST_BETA_UI, 2 , mean )
  # MEAN_OMEGA_UI =  apply( POST_OMEGA_UI, 2 , mean )
  Stage_Est_1 <- (Sys.time() - start)
  outputs <- list(POST_BETA_UI, Stage_Est_1)
  names(outputs) <- list("postbeta", "times")
  return(outputs)
}
