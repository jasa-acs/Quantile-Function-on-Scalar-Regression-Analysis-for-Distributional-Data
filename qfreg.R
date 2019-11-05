


####  qfreg.R is a code to fit the quantile functional regression model

### Input
### Argument		Type						Description
### X	N by A matrix	Covariates (N:# of observations)
### raw.dataset0		List of length N				Raw dataset (each list includes the repeated measurement for the outcome).
### prebasis		P by K matrix				Pre-quantlets basis functions, which implies the selected beta cdf functions. (P: desired # of probability grids, and K: # of basis).
### quantlets		P by K matrix				Quantlets basis function to transform back to the original data space (P: desired # of probability grids, and K: # of basis).
### n.iter			Positive integer				# of MCMC iteration
### burn			Positive integer				# of burning
### ep1			Positive scalar				Prior for nu in the paper (inverse gamma dist)

### Output yields list of length 6  and say qfregfit  as the output object.
### Argument		Type						Description
### qfregfit[[1]]		M by KA matrix		MCMC samples for coefficients (A:desired # of covariates, K:desired # of basis functions, and M: # of MCMC iteration after burning.
### qfregfit[[2]]		A by K matrix		Estimations of intercept and covariates in the projected space.
### qfregfit[[3]]		Scalar			Running time.
### qfregfit[[4]]		K by N matrix		Empirical coefficients
### qfregfit[[5]]		P by K matrix		Basis function to transform back to the original data space (P: desired # of probability grids, and K: # of basis).
### qfregfit[[6]]		P by N matrix		Empirical quantile functions for each subject in data space.

qfreg <- function(X, Y.list, REDUCED_BASE = prebasis, orthobais = quantlets, n.iter = 2000, burn = 200, ep1 = 0.0064) {
  initial <- EmpQuant22(REDUCED_BASE, orthobais, X, Y.list) ## names(initial)

  mcmc <- MCMC_QUANTLET(X, Y = initial$sd_l2, r = initial$cluster, TB00 = initial$TB00, n.iter = n.iter, burn = burn, ep1 = ep1)

  postmean <- apply(mcmc[[1]], 2, mean)

  outputs <- list(
    mcmc[[1]], postmean, mcmc[[2]], initial$sd_l2,
    initial$sdPhi, initial$empQ
  )
  names(outputs) <- list("postsample", "postmean", "times", "empcoef", "phi", "Q")

  return(outputs)
}

####################################################################################

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


EmpQuant22 <- function(REDUCED_BASE, orthobais, X, Y.list, delta2 = 0.95, H = 7) {
  EmpCoefs_quant_fit <- EmpCoefs_reg(orthobais, Y.list)
  empCoef0 <- EmpCoefs_quant_fit$emp
  Q <- EmpCoefs_quant_fit$Q


  library(pracma)

  if (dim(REDUCED_BASE)[1] != dim(orthobais)[1]) {
    biggrids.long <- dim(orthobais)[1]
    set.seed(123)
    biggrids <- sort(sample(dim(REDUCED_BASE)[1], biggrids.long))
    REDUCED_BASE_4096 <- REDUCED_BASE[  biggrids, ]
  }
  if (dim(REDUCED_BASE)[1] == dim(orthobais)[1]) {
    REDUCED_BASE_4096 <- REDUCED_BASE
  }



  Gram <- gramSchmidt(REDUCED_BASE, tol = .Machine$double.eps^0.5)
  Gram$Q

  r <- Cluster(REDUCED_BASE, Gram$Q, H = 7)


  d_l2 <- empCoef0
  K <- dim(d_l2)[2]
  Tp <- dim(orthobais)[1]
  Phi <- EmpCoefs_quant_fit$Psi


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

  outputs <- list(sd_l2, sdPhi, scale_f, d_l2, Phi, TB00, r, empCoef0, Q)
  names(outputs) <- list(
    "sd_l2", "sdPhi", "scale_f", "d_l2", "Phi", "TB00", "cluster",
    "empcoef", "empQ"
  )
  return(outputs)
}
#########################################################################################################

####  EmpCoefs_reg  computes scaled empirical coefficients based on the list outcome

EmpCoefs_reg <- function(B, Y.list) {
  n <- length(Y.list)
  ngrid <- dim(B)[1]
  gridp <- seq(1 / (ngrid + 1), ngrid / (ngrid + 1), 1 / (ngrid + 1))

  Psi <- cbind(rep(1, ngrid), B)

  emp <- matrix(NA, nrow = n, ncol = (dim(B)[2] + 1))
  Qregular <- matrix(NA, nrow = ngrid, ncol = n)

  for (i in 1:n) { ### i = 1
    y <- Y.list[[i]]
    y.long <- length(y)
    grideach <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
    if (length(grideach) != ngrid) {
      Qregular[, i] <- qhat(y, py = grideach, pred.p = gridp)
    }
    if (length(grideach) == ngrid) {
      Qregular[, i] <- y
    }
    emp[i, ] <- solve(t(Psi) %*% Psi) %*% t(Psi) %*% Qregular[, i]
  }
  output <- list(emp, Psi, Qregular)
  names(output) <- list("emp", "Psi", "Q")
  return(output)
}
#########################################################################################################

####  qhat is a local function to compute subject-specific Q_i-hat on common discrite grid, where it works within qfreg.R

qhat <- function(y = Q[select ], py = pi, pred.p = P) {
  y <- sort(y)
  ## py = round(  seq( 1/( length( y)+1 ),  length(  y)/( length( y)+1 ) , 1/( length( y)+1 ) )  , 4)
  indice <- matrix(NA, ncol = 7, nrow = length(pred.p))
  var.py <- seq(length(py))
  for (i in 1:length(pred.p)) { ## i = 119
    indice[i, 1] <- ifelse(sum(py <= pred.p[i]) == 0, NA, max(var.py[ py <= pred.p[i] ]))
    indice[i, 2] <- ifelse(sum(py >= pred.p[i]) == 0, NA, min(var.py[ py >= pred.p[i] ]))
    if ((is.na(indice[i, 1]) == TRUE) | (is.na(indice[i, 2]) == TRUE)) {
      {
        indice[i, 3] <- indice[i, 1]
        indice[i, 4] <- indice[i, 2]
        indice[i, 5] <- y[ indice[i, 1] ]
        indice[i, 6] <- y[ indice[i, 2] ]
        indice[i, 7] <- mean(c(indice[i, 5], indice[i, 6]), na.rm = TRUE)
      }
    }
    if ((is.na(indice[i, 1]) == FALSE) & (is.na(indice[i, 2]) == FALSE)) {
      if (indice[i, 1] == indice[i, 2]) {
        indice[i, 7] <- y[indice[i, 1]]
      }
      if (indice[i, 1] != indice[i, 2]) {
        indice[i, 7] <- y[indice[i, 1]]
        indice[i, 3] <- (py[ indice[i, 2]] - pred.p[i]) / (py[ indice[i, 2] ] - py[ indice[i, 1] ])
        indice[i, 4] <- (pred.p[i] - py[ indice[i, 1] ]) / (py[ indice[i, 2] ] - py[ indice[i, 1] ])
        indice[i, 5] <- y[ indice[i, 1] ]
        indice[i, 6] <- y[ indice[i, 2] ]
        indice[i, 7] <- indice[i, 3] * indice[i, 5] + indice[i, 4] * indice[i, 6]
      }
      if ((py[ indice[i, 2]] == py[ indice[i, 1] ])) {
        indice[i, 7] <- max(y[ indice[i, 1] ], y[ indice[i, 2] ])
      }
    }
  }
  return(indice[, 7])
}

#########################################################################################################

####  Cluster is a local function to cluster basis indices based on their eigen-values within qfreg.R

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
#########################################################################################################

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
  g <- rep(100000, max(K, Px))

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

  POST_BETA_FULL <- MCMC_BETA[-c(1:burn), ]
  POST_OMEGA_FULL <- MCMC_OMEGA[-c(1:burn), ]
  POST_ALPHA <- MCMC_ALPHA[-c(1:burn), ]
  POST_TAU <- MCMC_TAU[-c(1:burn), ]

  outputs <- list(POST_BETA_FULL, Stage_Est_1)
  names(outputs) <- list("postbeta", "times")
  return(outputs)
}
###########################################################################################################
