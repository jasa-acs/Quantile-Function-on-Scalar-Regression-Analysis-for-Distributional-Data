


####  inference.R is the code to conduct the estimation and inference in data space based on the output produced from the qfreg function
#### Input
#### Argument	Type			Description
#### mcmcEmpCoef	M by KA matrix	MCMC samples for coefficients (A:desired # of covariates, K:desired # of basis functions, and
#### M: 		Scalar		# of MCMC iteration after burning.
#### BackTransfor	P by K matrix	Basis function to transform back to the original data space (P: desired # of probability grids, and K: # of basis).
#### X		N by A matrix	Covariates (N:# of observations)
#### signifit	Scalar in (0,1)	The level of size alpha (default=0.975 allow for the inference at the level of .05)
#### X1		N1 by A matrix	New covariates
#### p		Vector of length P  in (0,1)	Probability grids
#### n.sup		Positive integer 	# of grids points of density estimates
#### xranges	Vector of length 2 	Vector consisting of c(min, max) as the range of the domain for the density estimate.
#### Output
#### Argument	 Type			Description
#### est		A by K matrix	Estimations of intercept and covariates in the projected space.
#### estCIu		P by A matrix	Upper point confidence intervals for the intercept and covariates in the data space
#### estCIl		P by A matrix	Lower point confidence intervals for the intercept and covariates in the data space
#### DataEst	P by A matrix	Estimations of intercept and covariates in the data space.
#### jointCI	P by A by 2 array	JointCI[,,1] includes upper joint confidence intervals for the intercept and covariates in the data space whereas jointCI[,,2] includes lower joint confidence intervals for the intercept and covariates in the data space.
#### local_p	P by A-1 matrix	SimBas for covariates (logical values where true means significant and false means insignificant for the given level of size alpha).
#### global_p	Vector of length A-1	Posterior probability scores for covariates
#### mu_G		3 by N1 matrix	Point estimation of mu for new covariates X1 is given in the first row where the second and third rows give the upper and lower confidence intervals.
#### sigma_G	3 by N1 matrix	Point estimation of variance for new covariates X1 is given in the first row where the second and third rows give the upper and lower confidence intervals.
#### mu3_G		3 by N1 matrix	Point estimation of skewness for new covariates X1 is given in the first row where the second and third rows give the upper and lower confidence intervals.
#### mu4_G		3 by N1 matrix	Point estimation of kurtosis for new covariates X1 is given in the first row where the second and third rows give the upper and lower confidence intervals.
#### mu_diff	N1/2 by 10 matrix	Results for the mean different testing for two consecutive subjects (8th column means  posterior probability scores).
#### sigma_diff	N1/2 by 10 matrix	Results for the variance different testing for two consecutive subjects (8th column means  posterior probability scores).
#### mu3_diff	N1/2 by 10 matrix	Results for the  skewness  different testing for two consecutive subjects (8th column means  posterior probability scores).
#### mu4_diff	N1/2 by 10 matrix	Results for the  kurtosis  different testing for two consecutive subjects (8th column means  posterior probability scores).
#### den_G		n.sup-1 by N1 matrix	Density estimates for the new covariates.
#### runtime	Scalar		Running time.



inference <- function(mcmcEmpCoef, BackTransfor = NULL, X, signifit = 0.975,
                      X1 = NULL, p, n.sup = 77, xranges = c(-10, 60)) {
  BN <- dim(mcmcEmpCoef)[1]
  Px <- dim(X)[2]
  Tp <- dim(BackTransfor)[1]
  K <- dim(BackTransfor)[2]


  quantiles2 <- function(x, probs = c(1 - signifit, signifit)) {
    quantile(x, probs)
  }


  ## mcmcEmpCoef= mcmc_fit2[[1]]
  Sys.time() -> start
  ITF_BETA_FULL <- array(NA, c(BN, Tp, dim(X)[2])) ## dim( ITF_BETA_FULL)

  for (i in 1:(BN)) { ## i=1
    ITQ_BAYES <- BackTransfor %*% t(matrix(mcmcEmpCoef[ i, ], ncol = K, nrow = Px))
    for (j in 1:Px) {
      ITF_BETA_FULL[i, , j] <- ITQ_BAYES[, j]
    }
  }


  MEAN_BETA_FULL <- apply(mcmcEmpCoef, 2, mean)
  BETA_FULL <- matrix(MEAN_BETA_FULL, nrow = Px)

  QUAN_BETA_FULL <- apply(ITF_BETA_FULL, c(2, 3), quantiles2)
  QUAN_ucl_BETA_FULL <- QUAN_BETA_FULL[2, , ]
  QUAN_lcl_BETA_FULL <- QUAN_BETA_FULL[1, , ]

  DataSp_est <- BackTransfor %*% t(BETA_FULL)




  MCMC_F_BASIS_MEAN <- matrix(0, nrow = Tp, ncol = Px)
  MCMC_F_BASIS_VAR <- matrix(0, nrow = Tp, ncol = Px)
  MCMC_F_BASIS_CI <- array(0, c(Tp, Px, 2))
  MCMC_F_BASIS_Z <- matrix(0, nrow = BN, ncol = Px)


  for (j in 1:Px) { ## j =1
    M_j <- apply(ITF_BETA_FULL[, , j ], 2, mean)
    S_j <- (apply(ITF_BETA_FULL[, , j ], 2, var))^{
      0.5
    }
    Z_j <- (ITF_BETA_FULL[, , j ] - matrix(M_j, nrow = BN, ncol = Tp, byrow = TRUE)) %*% diag(S_j^{
      -1
    }, ncol = Tp, nrow = Tp)
    abs_Z_J <- apply(abs(Z_j), 1, max)
    qalpha <- quantile(abs_Z_J, signifit - 0.025)
    I_j_995 <- M_j + qalpha * S_j
    I_j_005 <- M_j - qalpha * S_j
    MCMC_F_BASIS_MEAN[, j ] <- M_j
    MCMC_F_BASIS_VAR[, j ] <- S_j
    MCMC_F_BASIS_Z[, j ] <- abs_Z_J ## msides

    MCMC_F_BASIS_CI[, j, 1 ] <- I_j_995
    MCMC_F_BASIS_CI[, j, 2 ] <- I_j_005
  }

  ##########################################################################################################
  PMAP_FULL0 <- matrix(0, nrow = Tp, ncol = Px - 1)
  for (j in 2:Px) { ## j =2
    M_j <- apply(ITF_BETA_FULL[, , j ], 2, mean)
    S_j <- (apply(ITF_BETA_FULL[, , j ], 2, var))^{
      0.5
    }
    for (i in 1:BN) {
      PMAP_FULL0[, (j - 1)] <- PMAP_FULL0[, (j - 1)] + (abs(M_j / S_j) <= MCMC_F_BASIS_Z[i, j]) + 0
    }
  }
  PMAP_FULL <- PMAP_FULL0 / BN

  flag <- function(x, probs = c((1 - signifit) * 2)) {
    (x <= (1 - signifit) * 2)
  }

  LPMAP_FULL <- apply(PMAP_FULL, 2, flag)
  GPMAP_FULL <- apply(PMAP_FULL, 2, min)
  ###############################################################################################################


  if (is.null(X1) != TRUE) {
    PX0 <- X1
    Px0 <- dim(X1)[2]
    ### BN=1800
    STAT_MU <- matrix(NA, nrow = BN, ncol = dim(PX0)[1])
    for (i in 1:BN) { ### i =1
      ITQ_BAYES0 <- matrix(NA, ncol = Px, nrow = Tp)
      for (j in 1:Px) {
        ITQ_BAYES0[, j  ] <- ITF_BETA_FULL[i, , j]
      }

      ITQ_BAYES <- ITQ_BAYES0 %*% t(PX0)
      STAT_MU[i, ] <- t(ITQ_BAYES) %*% rep(1 / Tp, Tp)
    }


    MU_BAYES <- rbind(apply(STAT_MU, 2, mean), apply(STAT_MU, 2, quantiles2))


    STAT_VAR <- matrix(NA, nrow = BN, ncol = dim(PX0)[1])
    for (i in 1:BN) { ### i =1
      ITQ_BAYES0 <- matrix(NA, ncol = Px, nrow = Tp)
      for (j in 1:Px) {
        ITQ_BAYES0[, j  ] <- ITF_BETA_FULL[i, , j]
      }

      ITQ_BAYES <- ITQ_BAYES0 %*% t(PX0)
      CT_ITQ_BAYES <- ITQ_BAYES - matrix(rep(MU_BAYES[1, ], each = Tp, nrow = Tp), ncol = dim(PX0)[1])
      STAT_VAR[i, ] <- (t(CT_ITQ_BAYES^2) %*% rep(1 / Tp, Tp))^{
        0.5
      }
    }

    VAR_BAYES <- rbind(apply(STAT_VAR, 2, mean), apply(STAT_VAR, 2, quantiles2))


    STAT_MU3 <- matrix(NA, nrow = BN, ncol = dim(PX0)[1])
    STAT_MU4 <- matrix(NA, nrow = BN, ncol = dim(PX0)[1])
    for (i in 1:BN) { ### i =1
      ITQ_BAYES0 <- matrix(NA, ncol = Px, nrow = Tp)
      for (j in 1:Px) {
        ITQ_BAYES0[, j  ] <- ITF_BETA_FULL[i, , j]
      }
      ITQ_BAYES <- ITQ_BAYES0 %*% t(PX0)
      CT_ITQ_BAYES <- ITQ_BAYES - matrix(rep(MU_BAYES[1, ], each = Tp), nrow = Tp, ncol = dim(PX0)[1])
      ST_ITQ_BAYES <- (ITQ_BAYES - matrix(rep(MU_BAYES[1, ], each = Tp), nrow = Tp, ncol = dim(PX0)[1])) %*% diag(VAR_BAYES[1, ]^{
        -1
      }, ncol = dim(PX0)[1], nrow = dim(PX0)[1])
      STAT_MU3[i, ] <- t(ST_ITQ_BAYES^3) %*% rep(1 / Tp, Tp)
      STAT_MU4[i, ] <- t(CT_ITQ_BAYES^4) %*% rep(1 / Tp, Tp) / VAR_BAYES[1, ]^4
    }
    MU3_BAYES <- rbind(apply(STAT_MU3, 2, mean), apply(STAT_MU3, 2, quantiles2))
    MU4_BAYES <- rbind(apply(STAT_MU4, 2, mean), apply(STAT_MU4, 2, quantiles2))

    ###############################################################################################################

    spx <- dim(PX0)[1]
    ## dim(STAT_MU)
    STAT_MU_DIFF <- matrix(NA, nrow = spx / 2, ncol = 10)
    k.n <- 0
    for (k in 1:(spx / 2)) { ###    k =6
      i <- 2 * (k - 1) + 1
      j <- 2 * (k - 1) + 2
      k.n <- k.n + 1
      TEST_STATS <- STAT_MU[, j ] - STAT_MU[, i]

      M_j <- mean(TEST_STATS)
      S_j <- var(TEST_STATS)^{
        0.5
      }
      Z_j <- (TEST_STATS - M_j) / S_j
      Quan_Z_j_u <- quantile(Z_j, 0.975)
      Quan_Z_j_l <- quantile(Z_j, 0.025)

      STAT_MU_DIFF[k.n, 1 ] <- i
      STAT_MU_DIFF[k.n, 2 ] <- j
      ## STAT_MU_DIFF[k.n, 3, mets]= round( TRUE_VAR[j] -  TRUE_VAR[i] , 2)
      STAT_MU_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
      STAT_MU_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
      STAT_MU_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
      STAT_MU_DIFF[k.n, 9] <- ifelse(STAT_MU_DIFF[k.n, 3] == 0, 0, 1)
      STAT_MU_DIFF[k.n, 10 ] <- min(mean(TEST_STATS > 0), mean(TEST_STATS < 0))
    }

    STAT_VAR_DIFF <- matrix(NA, nrow = spx / 2, ncol = 10)
    k.n <- 0
    for (k in 1:(spx / 2)) { ###    k =6
      i <- 2 * (k - 1) + 1
      j <- 2 * (k - 1) + 2
      k.n <- k.n + 1
      TEST_STATS <- STAT_VAR[, j] - STAT_VAR[, i]

      M_j <- mean(TEST_STATS)
      S_j <- var(TEST_STATS)^{
        0.5
      }
      Z_j <- (TEST_STATS - M_j) / S_j
      Quan_Z_j_u <- quantile(Z_j, 0.975)
      Quan_Z_j_l <- quantile(Z_j, 0.025)

      STAT_VAR_DIFF[k.n, 1] <- i
      STAT_VAR_DIFF[k.n, 2] <- j
      ## STAT_VAR_DIFF[k.n, 3]= round( TRUE_VAR[j] -  TRUE_VAR[i] , 2)
      STAT_VAR_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
      STAT_VAR_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
      STAT_VAR_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
      STAT_VAR_DIFF[k.n, 9] <- ifelse(STAT_VAR_DIFF[k.n, 3] == 0, 0, 1)
      STAT_VAR_DIFF[k.n, 10 ] <- min(mean(TEST_STATS > 0), mean(TEST_STATS < 0))
    }


    STAT_MU3_DIFF <- matrix(NA, nrow = spx / 2, ncol = 10)
    k.n <- 0
    for (k in 1:(spx / 2)) { ###    k =6
      i <- 2 * (k - 1) + 1
      j <- 2 * (k - 1) + 2
      k.n <- k.n + 1
      TEST_STATS <- STAT_MU3[, j] - STAT_MU3[, i ]

      M_j <- mean(TEST_STATS)
      S_j <- var(TEST_STATS)^{
        0.5
      }
      Z_j <- (TEST_STATS - M_j) / S_j
      Quan_Z_j_u <- quantile(Z_j, 0.975)
      Quan_Z_j_l <- quantile(Z_j, 0.025)

      STAT_MU3_DIFF[k.n, 1] <- i
      STAT_MU3_DIFF[k.n, 2] <- j
      ## STAT_MU3_DIFF[k.n, 3]= round( TRUE_MU3[j] -  TRUE_MU3[i] , 2)
      STAT_MU3_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
      STAT_MU3_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
      STAT_MU3_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
      STAT_MU3_DIFF[k.n, 9] <- ifelse(STAT_MU3_DIFF[k.n, 3] == 0, 0, 1)
      STAT_MU3_DIFF[k.n, 10 ] <- min(mean(TEST_STATS > 0), mean(TEST_STATS < 0))
    }



    STAT_MU4_DIFF <- matrix(NA, nrow = spx / 2, ncol = 10)
    k.n <- 0
    for (k in 1:(spx / 2)) { ###    k =6
      i <- 2 * (k - 1) + 1
      j <- 2 * (k - 1) + 2
      k.n <- k.n + 1

      TEST_STATS <- STAT_MU4[, j] - STAT_MU4[, i]

      M_j <- mean(TEST_STATS)
      S_j <- var(TEST_STATS)^{
        0.5
      }
      Z_j <- (TEST_STATS - M_j) / S_j
      Quan_Z_j_u <- quantile(Z_j, 0.975)
      Quan_Z_j_l <- quantile(Z_j, 0.025)

      STAT_MU4_DIFF[k.n, 1] <- i
      STAT_MU4_DIFF[k.n, 2] <- j
      ## STAT_MU4_DIFF[k.n, 3]= round( TRUE_MU4[j] -  TRUE_MU4[i] , 2)
      STAT_MU4_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
      STAT_MU4_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
      STAT_MU4_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
      STAT_MU4_DIFF[k.n, 9] <- ifelse(STAT_MU4_DIFF[k.n, 3] == 0, 0, 1)
      STAT_MU4_DIFF[k.n, 10] <- min(mean(TEST_STATS > 0), mean(TEST_STATS < 0))
    }


    #########################################################################################################

    ## n.sup = 77
    xdomain <- seq(xranges[1], xranges[2], length.out = n.sup)
    I_DENSITY_FULL <- array(NA, c(BN, n.sup, dim(X1)[1]))
    DENSITY_FULL <- array(NA, c(BN, (n.sup - 1), dim(X1)[1]))


    for (i in 1:BN) { ### i =1
      ITQ_BAYES0 <- matrix(NA, ncol = Px, nrow = Tp)
      for (j in 1:Px) {
        ITQ_BAYES0[, j  ] <- ITF_BETA_FULL[i, , j]
      }
      ITQ_BAYES <- ITQ_BAYES0 %*% t(X1)

      for (h in 1:dim(X1)[1]) { ### h =1

        for (k in 1:n.sup) { ##  k=1
          if (sum(ITQ_BAYES[, h] <= xdomain[k]) == 0) {
            I_DENSITY_FULL[i, k, h] <- 0
          }
          if ((sum(ITQ_BAYES[, h] <= xdomain[k]) != 0)) {
            I_DENSITY_FULL[i, k, h] <- max(p[ ITQ_BAYES[, h] <= xdomain[k] ])
          }

          DENSITY_FULL[i, , h] <- diff(I_DENSITY_FULL[i, , h], 1) / diff(xdomain, 1)
        }
      } ### close
    } ## close where is.null(X1) ==FALSE


    MEAN_DENS_FULL <- apply(DENSITY_FULL, c(2, 3), mean)
    QUNT_DENS_FULL <- apply(DENSITY_FULL, c(2, 3), quantiles2)

    #########################################################################################
  } ## if closed!
  if (is.null(X1) == TRUE) {
    MU_BAYES <- NULL
    VAR_BAYES <- NULL
    MU3_BAYES <- NULL
    MU4_BAYES <- NULL
    STAT_MU_DIFF <- NULL
    STAT_VAR_DIFF <- NULL
    STAT_MU3_DIFF <- NULL
    STAT_MU4_DIFF <- NULL
    QUNT_DENS_FULL <- NULL
    MEAN_DENS_FULL <- NULL
  }

  Stage_Est_1 <- (Sys.time() - start)

  outputs <- list(
    BETA_FULL, QUAN_ucl_BETA_FULL, QUAN_lcl_BETA_FULL, DataSp_est,
    MCMC_F_BASIS_CI, LPMAP_FULL, GPMAP_FULL, MU_BAYES, VAR_BAYES, MU3_BAYES, MU4_BAYES,
    STAT_MU_DIFF, STAT_VAR_DIFF, STAT_MU3_DIFF, STAT_MU4_DIFF, QUNT_DENS_FULL, MEAN_DENS_FULL, Stage_Est_1
  )

  names(outputs) <- list(
    "est", "estCIu", "estCIl", "DataEst",
    "jointCI", "local_p", "global_p", "mu_G", "sigma_G", "mu3_G", "mu4_G",
    "mu_diff", "sigma_diff", "mu3_diff", "mu4_diff", "denCI_G", "den_G", "runtime"
  )

  return(outputs)
}

#########################################################################################################


Cluster <- function(REDUCED_BASE, orthobais, H = 7) {
  lambda_v <- t(REDUCED_BASE) %*% orthobais
  est_eigen <- round(diag(t(lambda_v) %*% lambda_v), 8)
  if (length(est_eigen) >= 6) {
    H1 <- H - 1
    hc <- hclust(dist(est_eigen)^2, "cen")
    r <- cutree(hc, k = H1)
  }
  if (length(est_eigen) < 6) {
    H1 <- H - 1
    r <- seq(length(est_eigen))
  }
  return(c(1, r + 1))
}



## r    = cluster indicator ( vector)
## B00  = initial leatst square est (matrix)
## A00  = zeros- non zero indicator for B00 (matrix)
## X    = covariates ( N x Px )
## Y  = scaled empirical coefficients( n x K)






####  MCMC_QUANTLET is a local function to run MCMC within qfreg.R

##  MCMC_QUANTLET-- conducts MCMC computaions
##  input
## 	X: covariates
## 	Y: empirical coefficients
## 	r: clustering indicators
## 	TB00: initial values of functional coefficients
## 	n.iter: # of MCMC iteration
## 	burn:   # of burning
## 	ep1:    prior for nu in the paper (inverse gamma dist)
##  output
## 	list[[1]]: MCMC samples
## 	list[[2]]: running time


##  EmpCoefs -- computes empirical coefficients
##  input
## 	B: basis functions
## 	Y.list: raw data list
##  output
## 	matrix consisting of empirical coefficients

MCMC_QUANTLET <- function(X, Y = sd_l2, r, TB00, n.iter = 2000, burn = 200, zeroKept = FALSE, ep1 = 0.0064) {

  ###  set up  ###
  N <- dim(X)[1]
  Px <- dim(X)[2]
  K <- dim(Y)[2]

  H <- length(unique(r))
  B00 <- solve(t(X) %*% X) %*% t(X) %*% Y
  A00 <- TB00
  A00[ A00 == 0  ] <- 0
  A00[ A00 != 0  ] <- 1
  C00 <- matrix(rep(r, time = Px), nrow = Px, ncol = K, byrow = TRUE)
  vec_C00 <- as.vector(C00)
  PA00 <- matrix(NA, ncol = K, nrow = Px)
  vec_PA00 <- as.vector(PA00)

  SSE_K <- diag(t(Y - X %*% B00) %*% (Y - X %*% B00))
  sigma2_K <- SSE_K / N

  nu0 <- ep1
  g <- rep(100000, K)

  one <- 1
  ## Y=sd_l2

  #################################### Fully Bayesian#################################################################
  Sys.time() -> start

  ## n.iter    ##000

  set.seed(n.iter + 1981)

  MCMC_BETA <- matrix(0, ncol = Px * K, nrow = n.iter)
  MCMC_OMEGA <- matrix(0, ncol = K, nrow = n.iter)

  MCMC_ALPHA <- matrix(0, ncol = Px * K, nrow = n.iter)


  MCMC_PHI <- matrix(NA, ncol = length(unique(r)), nrow = burn)
  MCMC_GAM <- matrix(NA, ncol = length(unique(r)), nrow = burn)
  ## MCMC_PHI  = matrix(  NA, ncol= H ,  nrow =  n.iter  )
  ## MCMC_GAM  = matrix(  NA, ncol= H ,  nrow =  n.iter  )

  ## MCMC_OB  = matrix(  0, ncol= Px*K  ,  nrow =  n.iter  )
  MCMC_PA <- matrix(0, ncol = Px * K, nrow = burn)
  MCMC_TAU <- matrix(0, ncol = Px * K, nrow = burn)
  PA00 <- matrix(NA, ncol = K, nrow = Px)
  TAU00 <- matrix(NA, ncol = K, nrow = Px)




  for (it in 1:n.iter) { ###   it =1


    if (it <= burn) {
      ## new_x = rbinom(N, 1, Ft  )
      ## X[, 13] =new_x

      if (it == 1) {
        sigma2_K0 <- sigma2_K
        temp.B00 <- B00
        MCMC_BETA[ 1, ] <- as.vector(B00)
        temp.A00 <- A00
        temp.C00 <- C00
        MCMC_PA[ 1, ] <- as.vector(PA00)
        temp.PA00 <- PA00
        temp.TAU00 <- TAU00
      }
      if (it != 1) {
        sigma2_K0 <- MCMC_OMEGA[(it - 1), ]
        temp.B00 <- matrix(MCMC_BETA[(it - 1), ], ncol = K, nrow = Px)
        temp.A00 <- matrix(MCMC_ALPHA[(it - 1), ], ncol = K, nrow = Px)
      }



      ############## Empirical Bayes part  ##############

      for (j in 1:Px) { ## j =4;   j =1 ;  j =2 ;   j =3 ;

        if ((Px - 1) != 1) {
          B_j_ <- matrix(temp.B00[ -j, ], ncol = K, nrow = (Px - 1))
          X_j_ <- X[, -j]
        }
        if ((Px - 1) == 1) {
          B_j_ <- matrix(temp.B00[ -j, ], ncol = K, nrow = (Px - 1), byrow = TRUE)
          X_j <- matrix(X[, -j], ncol = (Px - 1), nrow = N, byrow = TRUE)
        }

        B_j <- temp.B00[  j, ]
        X_j <- X[, j]

        V_j_all_h <- (sum(X_j^2) / sigma2_K0)^{
          -1
        }
        eta_j_all_h <- B_j / sqrt(V_j_all_h)

        # cluster by engen  ##indicator by trashold by energy!
        cluster_j <- C00[j, ] ##   cluster_j_all_h
        nonzeor_j <- A00[j, ] ##   indicator_j_all_h

        ##

        Ob_j_all_h <- rep(NA, K)

        for (h in 1:length(unique(r))) { ###  h =6 ## r
          hh <- unique(r)[h]

          crruent_cluster_j <- (cluster_j == hh)

          nonzero_in_crruent_cluster_j <- (crruent_cluster_j) & (nonzeor_j == 1)
          zero_in_crruent_cluster_j <- (crruent_cluster_j) & (nonzeor_j == 0)

          total.obs_in_crruent_cluster_j <- sum(crruent_cluster_j)
          effective.obs_in_crruent_cluster_j <- sum(nonzero_in_crruent_cluster_j)
          prob_nonzero_in_crruent_cluster_j <- effective.obs_in_crruent_cluster_j / total.obs_in_crruent_cluster_j

          freedom_in_crruent_cluster_j <- effective.obs_in_crruent_cluster_j - 1

          MCMC_PHI[it, h] <- prob_nonzero_in_crruent_cluster_j


          if (prob_nonzero_in_crruent_cluster_j == 1) {
            nonzeor_j[ crruent_cluster_j] <- 1
            temp.TAU00[j, crruent_cluster_j ] <- Inf
            MCMC_GAM[it, h] <- Inf
            temp.PA00[j, crruent_cluster_j ] <- 1
          }

          if (prob_nonzero_in_crruent_cluster_j == 0) {
            if (zeroKept == TRUE) {
              nonzeor_j[ crruent_cluster_j] <- 0
              temp.TAU00[j, crruent_cluster_j ] <- 0
              MCMC_GAM[it, h] <- 0
              temp.PA00[j, crruent_cluster_j ] <- 0
            }
            if (zeroKept == FALSE) {
              nonzeor_j[ crruent_cluster_j] <- 1
              temp.TAU00[j, crruent_cluster_j ] <- Inf
              MCMC_GAM[it, h] <- Inf
              temp.PA00[j, crruent_cluster_j ] <- 1
            }
          }



          if ((prob_nonzero_in_crruent_cluster_j != 0) & (prob_nonzero_in_crruent_cluster_j != 1)) {
            if (freedom_in_crruent_cluster_j >= 1) {
              G_j_h <- max(0, sum((eta_j_all_h^2)[  nonzero_in_crruent_cluster_j ]) / (sum(nonzero_in_crruent_cluster_j) - 1))
              MCMC_GAM[it, h] <- G_j_h

              temp.TAU00[ j, crruent_cluster_j  ] <- V_j_all_h[  crruent_cluster_j ] * G_j_h
              prob_j_h <- prob_nonzero_in_crruent_cluster_j

              if (sum(nonzero_in_crruent_cluster_j) >= 1) {
                prob0 <- prob_j_h / (1 - prob_j_h) / sqrt(1 + G_j_h) * exp(0.5 * eta_j_all_h[  nonzero_in_crruent_cluster_j ]^2 * (G_j_h / (1 + G_j_h)))
                Ob_j_all_h[  nonzero_in_crruent_cluster_j] <- ifelse(prob0 != Inf, prob0, 10^128)
                temp.PA00[j, nonzero_in_crruent_cluster_j] <- Ob_j_all_h[ nonzero_in_crruent_cluster_j] / (Ob_j_all_h[ nonzero_in_crruent_cluster_j] + 1)
              }
              if (sum(zero_in_crruent_cluster_j) >= 1) {
                Ob_j_all_h[   zero_in_crruent_cluster_j ] <- prob_j_h / (1 - prob_j_h) / sqrt(1 + G_j_h) * exp(0.5 * eta_j_all_h[  zero_in_crruent_cluster_j  ]^2 * (G_j_h / (1 + G_j_h)))
                temp.PA00[j, zero_in_crruent_cluster_j] <- Ob_j_all_h[ zero_in_crruent_cluster_j ] / (Ob_j_all_h[ zero_in_crruent_cluster_j] + 1)
              }
            } ## close  freedom_in_crruent_cluster_j > 1
            if (freedom_in_crruent_cluster_j == 0) {
              G_j_h <- Inf
              MCMC_GAM[it, h] <- G_j_h
              temp.TAU00[ j, crruent_cluster_j  ] <- V_j_all_h[  crruent_cluster_j ] * G_j_h
              prob_j_h <- prob_nonzero_in_crruent_cluster_j

              if (sum(nonzero_in_crruent_cluster_j) >= 1) {
                Ob_j_all_h[  nonzero_in_crruent_cluster_j] <- Inf
                temp.PA00[j, nonzero_in_crruent_cluster_j] <- 1
              }
              if (sum(zero_in_crruent_cluster_j) >= 1) {
                Ob_j_all_h[  zero_in_crruent_cluster_j] <- Inf
                temp.PA00[j, zero_in_crruent_cluster_j] <- 1
              }
            } ## close  freedom_in_crruent_cluster_j <= 1
          } ## close probability condition
        } ## close h
      } ## close j

      MCMC_TAU[it, ] <- as.vector(temp.TAU00)
      MCMC_PA [it, ] <- as.vector(temp.PA00)
    } ### close   if( it <= burn )



    ############## Main Estimation  ##############


    for (j in 1:Px) { ## j =2;   j =1 ;

      for (k in 1:K) { ## k =26



        if ((Px - 1) != 1) {
          B_j_ <- matrix(temp.B00[ -j, ], ncol = K, nrow = (Px - 1))
          X_j_ <- X[, -j]
        }
        if ((Px - 1) == 1) {
          B_j_ <- matrix(temp.B00[ -j, ], ncol = K, nrow = (Px - 1), byrow = TRUE)
          X_j <- matrix(X[, -j], ncol = (Px - 1), nrow = N, byrow = TRUE)
        }

        B_j <- temp.B00[  j, ]
        X_j <- X[, j]


        hat_xSx <- sum(X_j^2 / sigma2_K0[k])
        hat_xSy1 <- sum(X_j * Y[, k] / sigma2_K0[k])
        hat_xSy2 <- sum(X_j * (X_j_ %*% B_j_[, k]) / sigma2_K0[k]) #### t(X_j)%*%X_j_%*% B_j_[, k]/sigma2_K0[k]
        hat_xSy <- hat_xSy1 - hat_xSy2


        if (temp.TAU00[j, k] <= 0.000000001) {
          Vb <- 0
        }
        if (temp.TAU00[j, k] > 0.000000001) {
          Vb <- g[j] * solve(hat_xSx + 1 / temp.TAU00[j, k]) / (g[j] + 1)
        } ##

        Eb <- Vb %*% (hat_xSy)
        E <- matrix(rnorm(one, 0, 1), nrow = 1, ncol = one)

        if (temp.TAU00[j, k] <= 0.000000001) {
          beta_j <- t(c(Eb))
        }
        if (temp.TAU00[j, k] > 0.000000001) {
          beta_j <- t(t(E %*% chol(Vb)) + c(Eb))
        }

        u <- runif(1, 0, 1)
        if (u <= temp.PA00[j, k]) {
          beta_j <- beta_j
        }
        if (u > temp.PA00[j, k]) {
          beta_j <- 0
        }
        ### posit =   (k-1)*Px + j
        ### MCMC_BETA[ it,  posit ] = beta_j
        temp.B00[j, k ] <- beta_j
      } ##   close ## j =1
    } ##  close ## k =1
    MCMC_BETA[it, ] <- as.vector(temp.B00)
    temp.A00 <- temp.B00
    temp.A00 [  temp.A00 == 0  ] <- 0
    temp.A00 [  temp.A00 != 0  ] <- 1
    MCMC_ALPHA[it, ] <- as.vector(temp.A00)


    est <- matrix(MCMC_BETA[it, ], nrow = Px, ncol = K)
    SSE <- diag(t(Y - X %*% est) %*% (Y - X %*% est))
    MCMC_OMEGA[it, ] <- 1 / rgamma(K, (nu0 + N) / 2, (nu0 + SSE) / 2)
  }

  Stage_Est_1 <- (Sys.time() - start)
  ###########################################################################################################

  POST_BETA_FULL <- MCMC_BETA[-c(1:burn), ]
  POST_OMEGA_FULL <- MCMC_OMEGA[-c(1:burn), ]
  POST_ALPHA <- MCMC_ALPHA[-c(1:burn), ]
  POST_TAU <- MCMC_TAU[-c(1:burn), ]

  outputs <- list(POST_BETA_FULL, Stage_Est_1)
  names(outputs) <- list("postbeta", "times")
  return(outputs)
}
###########################################################################################################

##  MCMC_NRPCT  -- computes Gaussianity
##  input
## 	mcmcEmpCoef: MCMC samples for coefficients
## 	signifit: the level of size alpha (default=0.975 allow for the inference at the level of .05)
## 	X1: new covariates
##  output
## 	MEAN_NORMPC_X:  average of Gaussianity percentage for X1
## 	Q005_NORMPC_X:	lower confidence intervals for X1
## 	Q995_NORMPC_X:	upper confidence intervals for X1


MCMC_NRPCT <- function(mcmcEmpCoef, signifit = 0.975, X1) {
  Px <- dim(X1)[2]
  MCMC_PCT <- matrix(NA, nrow = dim(mcmcEmpCoef)[1], ncol = dim(X1)[1])
  for (i in 1:(dim(mcmcEmpCoef)[1])) { ## i=1
    basis.fit <- X1 %*% matrix(mcmcEmpCoef[ i, ], nrow = Px)
    basis.fit_2 <- basis.fit^2
    MCMC_PCT[i, ] <- apply(basis.fit_2[, 1:2], 1, sum) / apply(basis.fit_2, 1, sum)
  }
  quantiles2 <- function(x, probs = c(1 - signifit, signifit)) {
    quantile(x, probs)
  }

  MEAN_NORMPC_X <- round(apply(MCMC_PCT, 2, mean), 3)
  QUAN_NORMPC_X <- apply(MCMC_PCT, 2, quantiles2)
  Q995_NORMPC_X <- QUAN_NORMPC_X[2, ]
  Q005_NORMPC_X <- QUAN_NORMPC_X[1, ]

  normaltab <- round(cbind(MEAN_NORMPC_X, Q005_NORMPC_X, Q995_NORMPC_X), 3)
  return(normaltab)
}
#########################################################################################################

##  PREMONOTN -- computes epsilon-monotonicity
##  input
## 	fitQ:   predicted quantile functions
## 	cutoff: epsilon
##  output
## 	DIFF_FITS:  the first derivative on all regions p
## 	MONO_FITS:	the rate of epsilon-monotonicity
## 	stats:	basic statistics for the rate of epsilon-monotonicity for each predicted quantile


PREMONOTN <- function(fitQ, cutoff) {
  DIFF_FITS <- matrix(NA, nrow = (dim(fitQ)[1] - 1), ncol = dim(fitQ)[2])
  MONO_FITS <- rep(0, dim(fitQ)[2])
  for (i in 1:dim(fitQ)[2]) {
    DIFF_FITS[, i] <- diff(fitQ[, i], lag = 1)
    MONO_FITS[i] <- sum(DIFF_FITS[, i] >= 0) / (dim(fitQ)[1] - 1)
  }
  stats <- c(mean(MONO_FITS), median(MONO_FITS), min(MONO_FITS), sum(MONO_FITS >= cutoff))
  outputs <- list(DIFF_FITS, MONO_FITS, stats)
  return(outputs)
}

#########################################################################################################

##  NEW_X -- computes new covariates for a broad range of covariate values

NEW_X <- function(input.x, seednum = b) {
  set.seed(301 + b)
  new.X <- matrix(NA, ncol = dim(input.x)[2], nrow = dim(input.x)[1])
  for (i in 1:dim(input.x)[1]) { ## i=1
    temp.id <- sample(seq(dim(input.x)[1]), dim(input.x)[2], replace = TRUE)
    for (j in 1:dim(input.x)[2]) {
      new.X[i, j] <- input.x[ temp.id[j], j ]
    }
  }
  return(new.X)
}


#########################################################################################################

##  EmpCoefs -- computes empirical coefficients
##  input
## 	B: basis functions
## 	Y.list: raw data list
##  output
## 	matrix consisting of empirical coefficients

EmpCoefs <- function(B, Y.list = raw.dataset) {
  n <- length(Y.list)
  mats <- matrix(NA, nrow = n, ncol = (dim(B)[2] + 1))
  for (i in 1:n) { ### i = 1
    y <- Y.list[[i]]
    y.long <- length(y)
    set.seed(123 + i)
    grids <- sort(sample(dim(B)[1], y.long))
    Psi <- cbind(rep(1, length(grids)), B[ grids, ])
    mats[i, ] <- solve(t(Psi) %*% Psi) %*% t(Psi) %*% y
  }
  return(mats)
}


#########################################################################################################


##  EmpQuant2 -- computes scaled empirical coefficients, clustering, and initial coefficients
##  input
## 	empCoef0: empirical coefficients (it may be output of EmpCoefs function)
## 	REDUCED_BASE: reduced basis function
## 	orthobais: quantlets basis function
## 	X: covariates
## 	delta2: cutoff percentage of the energy for the functional coefficients
## 	H: # of clustering groups
##  output
## 	sd_l2: scaled quantlets empirical coefficients
## 	sdPhi: scaled quantlets basis function
## 	scale_f: scale factors
## 	d_l2: quantlets basis function
## 	Phi: quantlets basis function
## 	TB00: initial values of functional coefficients for the given delta2
## 	r: clustering indicators for the given H

EmpQuant2 <- function(empCoef0, REDUCED_BASE, orthobais, X, delta2 = 0.95, H = 7) {
  if (dim(REDUCED_BASE)[1] != dim(orthobais)[1]) {
    biggrids.long <- dim(orthobais)[1]
    set.seed(123)
    biggrids <- sort(sample(dim(REDUCED_BASE)[1], biggrids.long))
    REDUCED_BASE_4096 <- REDUCED_BASE[  biggrids, ]
  }
  if (dim(REDUCED_BASE)[1] == dim(orthobais)[1]) {
    REDUCED_BASE_4096 <- REDUCED_BASE
  }


  upper_R_T <- t(REDUCED_BASE_4096) %*% orthobais # dim(lambda_v)  ## each column lambda_v
  singluar <- round(diag(upper_R_T), 8) ## upper_R_T%*%t(upper_R_T)
  est_eigen <- c(200, singluar)

  if (length(est_eigen) >= 7) {
    H1 <- H
    hc <- hclust(dist(est_eigen)^2, "cen")
    r <- cutree(hc, k = H1)
  }

  if (length(est_eigen) < 7) {
    H1 <- H
    r <- seq(length(est_eigen))
  }


  d_l2 <- empCoef0
  K <- dim(d_l2)[2]
  Tp <- dim(orthobais)[1]

  Phi <- cbind(rep(1, Tp), orthobais)


  scale_f <- diag(t(d_l2) %*% d_l2)^{
    -0.5
  }
  sd_l2 <- d_l2 %*% diag(scale_f, K, K) ## diag( t(sd_l2)%*%sd_l2 )
  sdPhi <- Phi %*% diag(scale_f^{
    -1
  }, K, K)


  Px <- dim(X)[2]

  B00 <- solve(t(X) %*% X) %*% t(X) %*% sd_l2
  BETA_LSE <- B00
  ## 1->smallest!  64-> largest!

  ### ORD_BETA_LSE = matrix(NA, ncol= Px, nrow= K )   ## j=1

  ### for( j in 1:Px){  ORD_BETA_LSE[, j] = order( abs( BETA_LSE[j,] )  )  }

  JP_XI <- matrix(NA, nrow = Px, ncol = K)
  JP_XE <- matrix(NA, nrow = Px, ncol = K) ##  d_0_2[i,][JP_I[i,] ]
  JP_XB <- matrix(NA, nrow = Px, ncol = K)
  B00_2 <- B00^2
  for (i in 1:Px) { ##  i =1
    JP_XI[i, ] <- c(1, 2, (order(B00_2[i, -c(1, 2) ], decreasing = TRUE) + 2))
    ## JP_XI[i,] =   c(1,   (order(    B00_2[i, -c(1) ] , decreasing = TRUE ) + 1 )  )
    JP_XE[i, ] <- cumsum(B00_2[i, JP_XI[i, ]]) / sum(B00_2[i, JP_XI[i, ]])
    JP_XB[i, JP_XI[i, ]  ] <- JP_XE[i, ]
  }

  if (Px != length(delta2)) {
    cuts <- rep(delta2, Px)
  }
  if (Px == length(delta2)) {
    cuts <- delta2
  }
  # cuts = c( 0.667, 0.667,  0.675, 0.652, 0.667, 0.652, 0.65  )
  Set.off <- vector("list", Px)
  Set.on <- vector("list", Px)
  for (i in 1:Px) { ## i=3
    Set.off[[i]] <- seq(K)[ JP_XB[i, ] > cuts[i]   ]
    Set.on[[i]] <- seq(K)[ JP_XB[i, ] <= cuts[i]   ]
  }



  TB00 <- matrix(0, nrow = Px, ncol = K)
  Lambda_TB00 <- matrix(0, nrow = K, ncol = Px)
  for (i in 1:Px) { ## i=3
    if (is.null(Set.off[[i]]) == TRUE) {
      Lambda_TB00[, i] <- 0
    }
    if (is.null(Set.off[[i]]) != TRUE) {
      Lambda_TB00[  Set.off[[i]], i] <- 100000
      Lambda_TB00[ -Set.off[[i]], i] <- 0
    }
    TB00[i, ] <- ifelse((Lambda_TB00[, i] == 0), B00[i, ], 0)
  }

  outputs <- list(sd_l2, sdPhi, scale_f, d_l2, Phi, TB00, r)
  names(outputs) <- list("sd_l2", "sdPhi", "scale_f", "d_l2", "Phi", "TB00", "cluster")
  return(outputs)
}
#########################################################################################################
