#########################################################################################################

### getquantlets.R is a source code to drive the quantlets basis functions

###  getquantlets.R

###  Input
### Argument	       Type	 			Description
### raw.dataset	       List of length N		Raw dataset (each list includes the repeated measurement for the outcome).
### p	                   Vector of length P  in (0,1)	Probability grids.
### opt			Character 			It should be chosen either regular or irregular, which means the length of the probability grids for each subject.
### cutoff			Scalar in (0,1)		Near-lossless values.

### Output yields  list of length 6 and say computQuat as the output object.

### Argument	      Type				Description
### computQuat[[1]]	P by K matrix	      Pre-quantlets basis functions, which implies the selected beta cdf functions. (P: desired # of probability grids, and K: # of basis).
### computQuat[[2]]	List of length N	      Empirical quantiles for each subject.
### computQuat[[3]]	Vector of length Z	LCCC value for each choice of the basis function (Z: the length of the possible choices for basis).
### computQuat[[4]]	List of length Z   	Corresponding basis columns for each possible choice.
### computQuat[[5]]	Vector of length Z	# of the basis for each possible choice.
### computQuat[[6]]	Vector of length Z	Frequency for each possible choice.

#########################################################################################################

getquantlets <- function(raw.dataset, p = p1024, opt = "regular", cutoff = 0.990) {
  library(glmnet)
  library(MASS)


  decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
      nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }

  if (decimalplaces(log(length(p), base = 2)) != 0) {
    stop("length p must be 2 to the power")
  }
  n <- length(raw.dataset)

  if (opt == "irregular") {
    lasso.list <- vector("list", n)
    for (i in 1:n) { ## i =1
      y <- raw.dataset[[i]]
      y.long <- length(y)


      grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
      CDFBETA <- GENERATE_BETA_CDF(a1, a2, grid.p)
      NQ <- qnorm(grid.p, 0, 1)
      BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
      BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA)) 


      set.seed(12345 + i)
      lasso_fit <- glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE)
      cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3)

      zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
      selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros == FALSE ]
      lasso.list[[i]] <- selects
    }

    lasso.list1 <- lapply(lasso.list, CatNorm)
    lasso.nonzero.obs1 <- lapply(lasso.list1, rep.list)
    lasso.counts.fit <- CountBasis(lasso.list1, lasso.nonzero.obs1)

    lasso_IncidenceVec_i_ <- vector("list", n)
    for (i in 1:n) {
      lasso_fitIncd <- IncidenceVec(lasso.list1[-i], lasso.nonzero.obs1[-i])
      lasso_IncidenceVec_i_[[i]] <- lasso_fitIncd
    }


    lccc.lasso <- LCCC.FIT_1(
      leaveout.list = lasso_IncidenceVec_i_, remain.counts = lasso.counts.fit[[3]],
      remain.basis = lasso.counts.fit[[2]], Y.list = raw.dataset, maxim = length(p)
    )
    quantlet.set <- LCCC.PLOT_1(plength = length(p), Values = lccc.lasso[[1]], Y.list = raw.dataset, output1 = lasso.counts.fit[[1]], output2 = lasso.counts.fit[[2]], output3 = lasso.counts.fit[[3]], cutoff)


    CDFBETA <- GENERATE_BETA_CDF(a1, a2, p)
    NQ <- qnorm(p)
    BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
    BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
    REDUCED_BASE <- BETA_BASE_TOTAL_2
  } ## close if opt


  if (opt == "regular") {
    lasso.list <- vector("list", n)
    Y.list <- vector("list", n)
    grid.p <- p
    quantiles <- function(x, probs = grid.p) {
      quantile(x, probs, type = 7)
    }

    ## if(betype=="cdf"){
    CDFBETA <- GENERATE_BETA_CDF(a1, a2, grid.p)
    NQ <- qnorm(grid.p, 0, 1)
    BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
    BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
    ##   } ## close cdf



    for (i in 1:n) { ## i =1
      y <- quantiles(raw.dataset[[i]])
      Y.list[[i]] <- y
      y.long <- length(y)
      set.seed(12345 + i)
      lasso_fit <- glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE)
      cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3)
      zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
      selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros == FALSE ]
      lasso.list[[i]] <- selects
    } ## close for i in 1:n


    lasso.list1 <- lapply(lasso.list, CatNorm)
    lasso.nonzero.obs1 <- lapply(lasso.list1, rep.list)
    lasso.counts.fit <- CountBasis(lasso.list1, lasso.nonzero.obs1)


    lasso_IncidenceVec_i_ <- vector("list", n)
    for (i in 1:n) {
      lasso_fitIncd <- IncidenceVec(lasso.list1[-i], lasso.nonzero.obs1[-i])
      lasso_IncidenceVec_i_[[i]] <- lasso_fitIncd
    }

    lccc.lasso <- LCCC.FIT(
      leaveout.list = lasso_IncidenceVec_i_, B = BETA_BASE_TOTAL_2, remain.counts = lasso.counts.fit[[3]],
      remain.basis = lasso.counts.fit[[2]], Y.list = Y.list, maxim = length(p)
    )

    quantlet.set <- LCCC.PLOT_1(plength = length(p), Values = lccc.lasso[[1]], Y.list = Y.list, output1 = lasso.counts.fit[[1]], output2 = lasso.counts.fit[[2]], output3 = lasso.counts.fit[[3]], cutoff)

    REDUCED_BASE <- BETA_BASE_TOTAL_2
  } ## close if opt



  plotinput <- list(REDUCED_BASE, Y.list, lccc.lasso[[1]], lasso.counts.fit[[1]], lasso.counts.fit[[2]], lasso.counts.fit[[3]])
  names(plotinput) <- list("REDUCED_BASE", "Y.list", "Values", "output1", "output2", "output3")

  return(plotinput)
}

##################################################################################################################
#########################################################################################################




### CatNorm is a source code to add 1 (Gaussian basis) to all individual basis set
### input  -> vector:basis set
### output -> vector:basis set including 1
CatNorm <- function(vector) {
  unique(c(0, 1, sort(vector)))
}

### rep.list is a source code to generate one vectors with the length of each individual basis set
### input  -> vector:basis set
### output -> vector:one vectors with the length of each individual basis set
rep.list <- function(vector) {
  rep(1, length(vector))
}


#########################################################################################################


### IncidenceVec is a source code to compute frequency for each selected basis
### input  -> setlistofobs :basis set(CatNorm),  frqlistofobs:frequency for each count(rep.list)
### output -> [[1]]:basis names across subjects, [[2]]:frequency for each basis across subjects

IncidenceVec <- function(setlistofobs, frqlistofobs) {
  colum <- unlist(setlistofobs) ## length(colum )
  obs <- unlist(frqlistofobs) ## length( obs)
  mytable_colum0 <- table(colum, obs)
  unique.colum <- as.numeric(rownames(margin.table(mytable_colum0, 1)))
  unique.count <- as.numeric(margin.table(mytable_colum0, 1))
  outputs <- list(unique.colum, unique.count)
  return(outputs)
}
#########################################################################################################


### CountBasis is a source code to generate nested subsets of basis, # basis in each nested subse, and  each count(# of votes) yeilding each nested basis
### input  -> nonzero.list: sparse set Di,       nonzero.obs: one vectors with the same length to Di
### output ->  [[1]]:the nested subsets of basis, [[2]]:# basis in each nested subset
###            [[3]]: each count(# of votes) yeilding each nested basis.

CountBasis <- function(nonzero.list, nonzero.obs) {
  colum <- unlist(nonzero.list) ## length(colum )
  obs <- unlist(nonzero.obs) ## length(obs)
  mytable_colum0 <- table(colum, obs)
  unique.colum <- as.numeric(rownames(margin.table(mytable_colum0, 1))) ## length( unique.colum )
  unique.count <- as.numeric(margin.table(mytable_colum0, 1)) ## length( unique.count )
  margin.count <- c(0, sort(unique(unique.count))) ## length(margin.count)
  long.set <- length(margin.count) - 1
  list.set <- vector("list", long.set)
  for (i in 1:long.set) {
    list.set[[i]] <- sort(unique.colum[ unique.count > margin.count[i] ])
  }
  outputs <- list(list.set, unlist(lapply(list.set, length)), margin.count[-length(margin.count)])
  return(outputs)
}
#########################################################################################################


## LCCC.FIT  --  computes the leave-one-out concordance correlation index
##  input
## 	leaveout.list: leave-one-out list (it may be output of IncidenceVec function)
## 	remain.counts:  vector containing c value in paper  (it may be output[[2]] of CountBasis)
## 	remain.basis :  vector containing k(c) value in paper  (it may be output[[3]] of CountBasis)

## 	Y.list:  raw.dataset
## 	maxim:   scalar value (length of probability grids)
##  output
## 	list[[1]]: predicted values ( length(probability grids) times n times length(c)  array output ) based on leaveout.list


LCCC.FIT <- function(leaveout.list, B = BETA_BASE_TOTAL_2, remain.counts, remain.basis, Y.list, maxim = 1000) {
  n <- length(leaveout.list)
  active.set <- (remain.basis < maxim)
  if (tail(remain.counts[active.set], 1) == n - 1) {
    remain.counts[ length(remain.counts)] <- n - 2
  }
  feasible.long <- sum(active.set)
  max.long <- max(unlist(lapply(Y.list, length)))
  checks <- matrix(NA, nrow = n, ncol = feasible.long)
  Values <- array(NA, c(max.long, n, feasible.long))


  for (i in 1:n) { ### i = 35
    y <- Y.list[[i]]
    y.long <- length(y)
    set.seed(123 + i)
    grids <- sort(sample(dim(B)[1], y.long))
    Psi <- cbind(rep(1, length(grids)), B[ grids, ])

    for (j in 1:feasible.long) { ###  j = 20

      colum_i_ <- leaveout.list[[i]][[1]]
      obs_i_ <- leaveout.list[[i]][[2]] ## length(IncidenceVec_i_)

      SET1_i_ <- sort(colum_i_[  obs_i_ > remain.counts[active.set ][j]   ])
      smPsi_i_ <- Psi[, (SET1_i_) + 1 ] ## dim(smPsi_i_)

      Values[(1:y.long), i, j] <- try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% t(smPsi_i_) %*% y) ## dim(Qy)
      checks[i, j] <- length(SET1_i_)
    }
  }
  outputs <- list(Values, checks)
}
#########################################################################################################

## LCCC.FIT1  --  computes the leave-one-out concordance correlation index
##  input
## 	leaveout.list: leave-one-out list (it may be output of IncidenceVec function)
## 	remain.counts:  vector containing c value in paper  (it may be output[[2]] of CountBasis)
## 	remain.basis :  vector containing k(c) value in paper  (it may be output[[3]] of CountBasis)

## 	Y.list:  raw.dataset
## 	maxim:   scalar value (length of probability grids)
##  output
## 	list[[1]]: predicted values ( length(probability grids) times n times length(c)  array output ) based on leaveout.list


LCCC.FIT_1 <- function(leaveout.list, remain.counts, remain.basis, Y.list, maxim = 1000) {
  n <- length(leaveout.list)
  active.set <- (remain.basis < maxim)
  if (tail(remain.counts[active.set], 1) == n - 1) {
    remain.counts[ length(remain.counts)] <- n - 2
  }
  feasible.long <- sum(active.set)
  max.long <- max(unlist(lapply(Y.list, length)))
  checks <- matrix(NA, nrow = n, ncol = feasible.long)
  Values <- array(NA, c(max.long, n, feasible.long))


  for (i in 1:n) { ### i = 35
    y <- Y.list[[i]]
    y.long <- length(y)
    grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
    CDFBETA <- GENERATE_BETA_CDF(a1, a2, grid.p)
    NQ <- qnorm(grid.p, 0, 1)
    BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
    BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
    Psi <- cbind(rep(1, length(grid.p)), BETA_BASE_TOTAL_2)

    for (j in 1:feasible.long) { ###  j = 20

      colum_i_ <- leaveout.list[[i]][[1]]
      obs_i_ <- leaveout.list[[i]][[2]] ## length(IncidenceVec_i_)

      SET1_i_ <- sort(colum_i_[  obs_i_ > remain.counts[active.set ][j]   ])
      smPsi_i_ <- Psi[, (SET1_i_) + 1 ] ## dim(smPsi_i_)

      Values[(1:y.long), i, j] <- try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% t(smPsi_i_) %*% y) ## dim(Qy)
      checks[i, j] <- length(SET1_i_)
    }
  }
  outputs <- list(Values, checks)
}


#########################################################################################################

## LCCC.PLOT  --  creates lccc plot (rho0 vs k(c) )
##  input
## 	plength: length of probability grids
## 	Values:  predicted values (it may be output of LCCC.FIT )
## 	output1: list of vector including the indices of remained basis correponding to removing the basis with the low frequency K(c) in the paper  (it may be output[[1]] of CountBasis)
## 	output2: remain.counts   (it may be output[[2]] of CountBasis)
## 	output3: remain.basis    (it may be output[[3]] of CountBasis)
## 	cutoff: 1-epsilon
##  output
## 	list[[1]]: quantlets basis function
## 	list[[2]]: lccc value
## 	lccc plot


LCCC.PLOT <- function(plength, Values, Y.list, output1, output2, output3, cutoff) {
  n <- length(Y.list)
  Y.mat <- matrix(NA, nrow = max(unlist(lapply(Y.list, length))), ncol = length(Y.list))
  for (i in 1:n) {
    Y.mat[(1:length(Y.list[[i]])), i] <- Y.list[[i]]
  }


  lasso.long <- dim(Values)[3]
  lasso.values <- matrix(NA, nrow = n, ncol = lasso.long)

  for (j in 1:lasso.long) {
    lasso.values[, j] <- agree.cci(Y.mat, Values[, , j])
  }

  lasso.Chary_i_ <- apply(lasso.values, 2, mean)
  lasso.Chary1_i_ <- apply(lasso.values, 2, min)
  ## define names
  lasso.x <- output2[ (output2 < plength) ][1:lasso.long]
  lasso.x1 <- output3[ (output2 < plength) ][1:lasso.long]

  texts1 <- expression(paste(bar(rho), " vs K"))
  texts2 <- expression(paste(rho^0, " vs K"))

  cex.lab <- 0.9 ### -> Size of labels for both x and y axis!
  cex.axis <- 0.65 ### -> Size of coordinates for both x and y axis!
  cex.main <- 0.7 ### -> Size of main topic ( USELESS since we will cut title )
  cex.in <- 0.5 ### -> Size of all in boxplot (Ex: outlier )

  yaxis.at <- c(-2, seq(0.85, 1.0, by = 0.05), 1.5)
  yaxis.lab <- c("", "0.85", "0.90", "0.95", "1.00", "")
  ylim <- c(0.84, 1.0)

  numbasis <- sort(unique(c(lasso.x)), decreasing = TRUE)
  plotxaxis <- seq(length(numbasis))

  xaxis.at <- c(plotxaxis[1] - 15, plotxaxis, tail(plotxaxis, 1) + 15)
  xaxis.lab <- c("", numbasis, "")
  xaxis.lab1 <- c("", lasso.x1, "")


  par(mfrow = c(1, 2), mar = c(4.5, 2.5, 3, 2))

  plot(0, type = "o", xlab = "", ylab = "", cex.lab = 1.1, cex.axis = 0.7, axes = FALSE, ylim = ylim, xlim = c(0, length(numbasis) + 1))
  points(plotxaxis, lasso.Chary1_i_, col = "red", pch = 19)
  lines(plotxaxis, lasso.Chary1_i_, col = "red")

  axis(side = 2, at = yaxis.at, label = yaxis.lab, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
  axis(side = 2, at = yaxis.at, label = rep("", length(yaxis.at)), line = -0.5, , tick = TRUE, las = 1)

  axis(side = 1, at = xaxis.at, label = xaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis)
  axis(side = 1, at = xaxis.at, label = rep("", length(xaxis.at)), line = -0.5, tick = TRUE)
  axis(side = 1, at = xaxis.at, label = xaxis.lab, line = 0.5, tick = FALSE, cex.axis = cex.axis)
  mtext("C", side = 1, line = 0.5, at = 0.5, col = "blue", cex = 0.9)
  mtext(expression(K[C]), side = 1, line = 1.7, at = 0.5, col = "blue", cex = 0.9)
  title(texts2)

  plot(0, type = "o", xlab = "", ylab = "", cex.lab = 1.1, cex.axis = 0.7, axes = FALSE, ylim = ylim, xlim = c(0, length(numbasis) + 1))
  points(plotxaxis, lasso.Chary_i_, col = "blue", pch = 19)
  lines(plotxaxis, lasso.Chary_i_, col = "blue")

  axis(side = 2, at = yaxis.at, label = yaxis.lab, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
  axis(side = 2, at = yaxis.at, label = rep("", length(yaxis.at)), line = -0.5, , tick = TRUE, las = 1)
  axis(side = 1, at = xaxis.at, label = xaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis)
  axis(side = 1, at = xaxis.at, label = rep("", length(xaxis.at)), line = -0.5, tick = TRUE)
  axis(side = 1, at = xaxis.at, label = xaxis.lab, line = 0.5, tick = FALSE, cex.axis = cex.axis)
  mtext("C", side = 1, line = 0.5, at = 0.5, col = "blue", cex = 0.9)
  mtext(expression(K[C]), side = 1, line = 1.7, at = 0.5, col = "blue", cex = 0.9)
  title(texts1)

  id1 <- lasso.Chary1_i_ >= cutoff
  quantlet.set <- output1[[  max(seq(sum(id1))) ]][-1]
  outputs <- list(quantlet.set, lasso.Chary1_i_[ max(seq(sum(id1))) ])
  names(outputs) <- list("quantlet", "lccc")
  return(outputs)
}
#########################################################################################################
## LCCC.PLOT1  --  creates lccc plot (rho0 vs k(c) )
##  input
## 	plength: length of probability grids
## 	Values:  predicted values (it may be output of LCCC.FIT )
## 	output1: list of vector including the indices of remained basis correponding to removing the basis with the low frequency K(c) in the paper  (it may be output[[1]] of CountBasis)
## 	output2: remain.counts   (it may be output[[2]] of CountBasis)
## 	output3: remain.basis    (it may be output[[3]] of CountBasis)
## 	cutoff: 1-epsilon
##  output
## 	list[[1]]: quantlets basis function
## 	list[[2]]: lccc value
## 	lccc plot

LCCC.PLOT_1 <- function(plength, Values, Y.list, output1, output2, output3, cutoff) {
  n <- length(Y.list)
  Y.mat <- matrix(NA, nrow = max(unlist(lapply(Y.list, length))), ncol = length(Y.list))
  for (i in 1:n) {
    Y.mat[(1:length(Y.list[[i]])), i] <- Y.list[[i]]
  }


  lasso.long <- dim(Values)[3]
  lasso.values <- matrix(NA, nrow = n, ncol = lasso.long)

  for (j in 1:lasso.long) {
    lasso.values[, j] <- agree.cci(Y.mat, Values[, , j])
  }

  lasso.Chary_i_ <- apply(lasso.values, 2, mean)
  lasso.Chary1_i_ <- apply(lasso.values, 2, min)
  ## define names
  lasso.x <- output2[ (output2 < plength) ][1:lasso.long]
  lasso.x1 <- output3[ (output2 < plength) ][1:lasso.long]

  texts1 <- expression(paste(bar(rho), " vs K"))
  texts2 <- expression(paste(rho^0, " vs K"))

  cex.lab <- 0.9 ### -> Size of labels for both x and y axis!
  cex.axis <- 0.65 ### -> Size of coordinates for both x and y axis!
  cex.main <- 0.7 ### -> Size of main topic ( USELESS since we will cut title )
  cex.in <- 0.5 ### -> Size of all in boxplot (Ex: outlier )

  yaxis.at <- c(-2, seq(0.85, 1.0, by = 0.05), 1.5)
  yaxis.lab <- c("", "0.85", "0.90", "0.95", "1.00", "")
  ylim <- c(0.84, 1.0)

  numbasis <- sort(unique(c(lasso.x)), decreasing = TRUE)
  plotxaxis <- seq(length(numbasis))

  xaxis.at <- c(plotxaxis[1] - 15, plotxaxis, tail(plotxaxis, 1) + 15)
  xaxis.lab <- c("", numbasis, "")
  xaxis.lab1 <- c("", lasso.x1, "")


  ## par( mfrow=c(1,1), mar=c(4.5, 2.5, 3,2 ))

  plot(0, type = "o", xlab = "", ylab = "", cex.lab = 1.1, cex.axis = 0.7, axes = FALSE, ylim = ylim, xlim = c(0, length(numbasis) + 1))
  points(plotxaxis, lasso.Chary1_i_, col = "red", pch = 19)
  lines(plotxaxis, lasso.Chary1_i_, col = "red")

  axis(side = 2, at = yaxis.at, label = yaxis.lab, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
  axis(side = 2, at = yaxis.at, label = rep("", length(yaxis.at)), line = -0.5, , tick = TRUE, las = 1)

  axis(side = 1, at = xaxis.at, label = xaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis)
  axis(side = 1, at = xaxis.at, label = rep("", length(xaxis.at)), line = -0.5, tick = TRUE)
  axis(side = 1, at = xaxis.at, label = xaxis.lab, line = 0.5, tick = FALSE, cex.axis = cex.axis)
  mtext("C", side = 1, line = 0.5, at = 0.5, col = "blue", cex = 0.9)
  mtext(expression(K[C]), side = 1, line = 1.7, at = 0.5, col = "blue", cex = 0.9)
  ## title(texts2 )

  id2 <- (lasso.x) <= (length(Y.list) - 1) ## length(lasso.x)
  id3 <- lasso.Chary1_i_ >= cutoff
  id1 <- id3 & id2
  if (sum(id1) == 0) {
    id <- id2
    this <- min(seq(length(id))[id])
  }
  if (sum(id1) != 0) {
    id <- id1
    this <- max(seq(length(id))[id])
  }

  quantlet.set <- output1[[  this   ]] [-1]
  outputs <- list(quantlet.set, lasso.Chary1_i_[ this ])
  names(outputs) <- list("quantlet", "lccc")
  return(outputs)
}

#########################################################################################################
## GENERATE_BETA_CDF --  creates beta basis functions
##  input
## 	alpha: vector containing sequence of beta parameter
## 	beta:  vector containing sequence of beta parameter
## 	index.p: probability grids on (0,1)
##  output
##  	matrix containing # of beta parameters times the length of index.p

GENERATE_BETA_CDF <- function(alpha, beta, index.p) {
  n1 <- as.character(alpha)
  n2 <- as.character(beta)
  BETASCDF0 <- matrix(NA, ncol = length(index.p), nrow = length(n2) * length(n1))
  for (i in 1:length(n1)) { ##   i=1;   j=12;   a1[i]    a2[j]
    for (j in 1:length(n2)) {
      rowth <- (j - 1) * length(n2) + i
      BETASCDF0[ rowth, ] <- pbeta(index.p, alpha[i], beta[j])
      BETASCDF0[ rowth, ] <- round(centering.function(BETASCDF0[ rowth, ], scale = TRUE), 7)
    }
  }
  ## name.mat = outer( paste0( "(", n1 , se="") ,paste0(",", n2, ")" , se="")  ,  FUN=paste ,sep="")
  ## matrix.rownames = as.vector(  name.mat )
  ## outputs= list( BETASCDF0 ,  name.mat , matrix.rownames )
  outputs <- BETASCDF0
  return(outputs)
}

#########################################################################################################

a1 <- c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1))
a2 <- c(seq(0.1, 1, by = 0.1), seq(2, 100, by = 1))

#########################################################################################################

## centering.function --  computes centering and scaling column by column
##  input
## 	raw.x: any matrix with n times p
## 	scale: FALSE or TRUE (doing only centering if scale=FALSE )
##  output: the matrix centerted or/and scaled

centering.function <- function(raw.x, scale = FALSE) {
  object.x <- as.matrix(raw.x)
  x.n <- dim(object.x)[1]
  one <- rep(1, x.n)
  meanx <- drop(one %*% object.x) / x.n
  mean1 <- as.matrix(meanx)
  n1 <- dim(mean1)[1]
  jj <- rep(1, x.n) %*% t(rep(1, dim(object.x)[2]))
  mean.mat <- diag(meanx, nrow = n1, ncol = n1)
  cen.x <- object.x - jj %*% mean.mat
  normx <- sqrt(drop(one %*% (cen.x^2)))
  normx1 <- as.matrix(normx)
  s1 <- dim(normx1)[1]
  scale.mat <- diag(1 / normx, nrow = s1, ncol = s1)
  cen.scale.x <- cen.x %*% scale.mat
  if (scale == FALSE) {
    return(cen.x)
  }
  if (scale == TRUE) {
    return(cen.scale.x)
  }
}
#########################################################################################################

## WaveletDenose-- conducts wavelet denoising procedure
##  input: any matrix with n times p
##  output: array[,,1] includes denoised values based on the usual wavelet method
## 	    array[,,2] includes denoised values based on nondecimal method

WaveletDenose <- function(object, filter.number = 2, family = "DaubExPhase") {
  object <- as.matrix(object)
  n <- dim(object)[1]
  p <- dim(object)[2]
  outputs <- array(NA, c(n, p, 2))
  for (i in 1:p) { ## i =1
    ywd <- wd(object[, i], filter.number = filter.number, family = family)
    FineCoefs <- accessD(ywd, lev = nlevelsWT(ywd) - 1)
    sigma <- mad(FineCoefs)
    utDJ <- sigma * sqrt(2 * log(1024))
    ywdT <- threshold(ywd, policy = "manual", value = utDJ)
    outputs[, i, 1] <- wr(ywdT)

    ywst <- wst(object[, i], filter.number = filter.number, family = family)
    FineWSTCoefs <- accessD(ywst, lev = nlevelsWT(ywst) - 1)
    sigmaWST <- mad(FineWSTCoefs)
    utWSTDJ <- sigmaWST * sqrt(2 * log(1024))
    ywstT <- threshold(ywst, policy = "manual", value = utWSTDJ)
    outputs[, i, 2] <- AvBasis(ywstT)
  }
  return(outputs)
}



##########################################################################
## concordance --  computes concordance correlation index
##  input
## 	yhat: fitted values
## 	y   : emprical true values
##  output: concordance value

concordance <- function(yhat, y) {
  2 * cov(yhat, y, use = "complete.obs") / (cov(yhat, yhat, use = "complete.obs") + cov(y, y, use = "complete.obs") + (mean(y, na.rm = TRUE) - mean(yhat, na.rm = TRUE))^2)
}
##########################################################################

## agree.cci--  computes concordance correlation index in the version of matrix
##  input
## 	input.yhat: fitted values as the matrix
## 	input.y   : emprical true values as the matrix
##  output: concordance value for the matrtix

agree.cci <- function(input.yhat, input.y) {
  k <- dim(input.yhat)[2]
  values <- rep(NA, k)
  for (i in 1:k) {
    values[i] <- concordance(input.y[, i], input.yhat[, i])
  }
  return(values)
}
##########################################################################

## uporder--  reorder basis columns based on the energy of its coefficients
##  input: d_0 is empirical coefficients
##  output: reordered basis columns in decreasing order


uporder <- function(d_0) {
  d_0_2 <- d_0^2
  tsum <- apply(d_0_2, 2, sum)
  reorder <- c(1, 2, (order(tsum [ -c(1, 2) ], decreasing = TRUE) + 2)) - 1
  return(reorder)
}


##########################################################################

## GENERATE_BETA_PDF --  creates beta basis functions
##  input
## 	alpha: vector containing sequence of beta parameter
## 	beta:  vector containing sequence of beta parameter
## 	index.p: probability grids on (0,1)
##  output
##  	matrix containing # of beta parameters times the length of index.p


GENERATE_BETA_PDF <- function(alpha = a1, beta = a2, index.p) {
  n1 <- as.character(alpha)
  n2 <- as.character(beta)
  BETASCDF0 <- matrix(NA, ncol = length(index.p), nrow = length(n2) * length(n1))
  for (i in 1:length(n1)) { ##   i=1;   j=12;   a1[i]    a2[j]
    for (j in 1:length(n2)) {
      rowth <- (j - 1) * length(n2) + i
      BETASCDF0[ rowth, ] <- dbeta(index.p, alpha[i], beta[j])
      BETASCDF0[ rowth, ] <- round(centering.function(BETASCDF0[ rowth, ], scale = TRUE), 7)
    }
  }
  ### name.mat = outer( paste0( "(", n1 , se="") ,paste0(",", n2, ")" , se="")  ,  FUN=paste ,sep="")
  ### matrix.rownames = as.vector(  name.mat )
  ### outputs= list( BETASCDF0 ,  name.mat , matrix.rownames )

  ########   dim( apply(  apply( BETASCDF0 , 1, is.na )) ## dim(apply( BETASCDF0 , 1, is.na ) )
  ########   idx = ( apply(  apply( BETASCDF0 , 1, is.na )    , 2 , sum )==0)
  outputs <- BETASCDF0[ -991, ] ## dim(BETASCDF0)
  return(outputs)
}
