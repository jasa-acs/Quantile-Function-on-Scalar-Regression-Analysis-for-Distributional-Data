

rm(list = ls())

library(glmnet)
library(MASS)
library(R.matlab)

library(pracma)
library(wavethresh)

source("getquantlets.R")
source("qfreg.R")
source("inference.R")
source("PcrQuant.R")




Sys.time() -> start
Q3cellfile <- readMat("Q3cell.mat")

raw.dataset <- vector("list", 64)
raw.length <- rep(0, 64)


for (i in 1:64) {
  raw.dataset[[i]] <- sort(as.numeric(unlist(Q3cellfile[[1]][[i]])))
  raw.length[[i]] <- length(raw.dataset[[i]])
}




########################  lasso selection! 30min########################

lasso.list <- vector("list", 64)
for (i in 1:64) {
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

Times.over.lasso <- (Sys.time() - start)


######################## Compute common basis  ########################

Sys.time() -> start
lasso.list1 <- lapply(lasso.list, CatNorm) ## 1 term fix
lasso.nonzero.obs1 <- lapply(lasso.list1, rep.list) ## D_i, generate 1 vector
lasso.counts.fit <- CountBasis(lasso.list1, lasso.nonzero.obs1)


n <- 64

lasso_IncidenceVec_i_ <- vector("list", n) ### generate leave-one-out
for (i in 1:n) {
  lasso_fitIncd <- IncidenceVec(lasso.list1[-i], lasso.nonzero.obs1[-i])
  lasso_IncidenceVec_i_[[i]] <- lasso_fitIncd
}



p1024 <- signif(seq(0.001, 0.999, length = 1024), 4)
  quantiles_p <- function(x, probs = p1024) {
    quantile(x, probs, type=6)
}


Qy = matrix( round(unlist( lapply( raw.dataset,  quantiles_p )  ),3) , 1024 )
write.csv(Qy, "Qy_1024.csv", row.names = FALSE )
Qy_3 <- read.csv("Qy_1024.csv", header = TRUE)
Qy <- Qy_3






leaveout.list <- lasso_IncidenceVec_i_
remain.counts <- lasso.counts.fit[[3]]
remain.basis <- lasso.counts.fit[[2]]
Y.list <- raw.dataset

active.set <- (remain.basis < 1000)
feasible.long <- sum(active.set)
n <- length(leaveout.list)

max.long <- max(unlist(lapply(Y.list, length)))
checks <- matrix(NA, nrow = n, ncol = feasible.long)
Values <- array(NA, c(max.long, n, feasible.long))



for (i in 1:n) {
  y <- raw.dataset[[i]]
  y.long <- length(y)
  grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
  CDFBETA <- GENERATE_BETA_CDF(a1, a2, grid.p)
  NQ <- qnorm(grid.p, 0, 1)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
  BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
  Psi <- cbind(rep(1, length(grid.p)), BETA_BASE_TOTAL_2)

  for (j in 1:feasible.long) {
    colum_i_ <- leaveout.list[[i]][[1]]
    obs_i_ <- leaveout.list[[i]][[2]]

    SET1_i_ <- sort(colum_i_[  obs_i_ > remain.counts[active.set ][j]   ])
    smPsi_i_ <- Psi[, (SET1_i_) + 1 ]

    Values[(1:y.long), i, j] <- try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% t(smPsi_i_) %*% y)
    checks[i, j] <- length(SET1_i_)
  }
}


##################### leave one-out cci ! for beta basis ##################################################################################
Y.mat <- matrix(NA, nrow = max.long, ncol = n)
for (i in 1:64) {
  Y.mat[(1:length(raw.dataset[[i]])), i] <- raw.dataset[[i]]
}


lasso.long <- sum(lasso.counts.fit[[2]] < 1000)
lasso.values <- matrix(NA, nrow = n, ncol = lasso.long)

for (j in 1:lasso.long) {
  lasso.values[, j] <- agree.cci(Y.mat, Values[, , j])
}

lasso.Chary_i_ <- apply(lasso.values, 2, mean)
lasso.Chary1_i_ <- apply(lasso.values, 2, min)

lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x1 <- lasso.counts.fit[[3]][ (lasso.counts.fit[[2]] < 1000) ]


##################### leave one-out cci ! for PCA ##################################################################################

## LCCC.FIT  --  computes the leave-one-out concordance correlation index based on PC regression
## TruncPCA  --  computes PCs basis with its maximum number (max.K)

fit.pca <- TruncPCA(Qy, length.train = NULL, max.K = 64)
lccc.pca <- LCCC.PCR.FIT(Y = Qy, B = fit.pca$v, remain.basis = lasso.counts.fit[[2]])



CDFBETA <- GENERATE_BETA_CDF(a1, a2, p1024)
NQ <- qnorm(p1024, 0, 1)
BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))


raw.dataset1 <- vector("list", 64)
for (i in 1:64) {
  raw.dataset1[[i]] <- Qy[, i]
}


##################### leave one-out cci ! for quanvelts ##################################################################################


lcccwave.lasso <- LCCC.WAVE.FIT(
  global.list = lasso.counts.fit[[1]], B = BETA_BASE_TOTAL_2,
  remain.counts = lasso.counts.fit[[3]],
  remain.basis = lasso.counts.fit[[2]], Y.list = raw.dataset1
)

Y.mat <- matrix(NA, nrow = 1024, ncol = n)
for (i in 1:64) {
  Y.mat[, i] <- raw.dataset1[[i]]
}


lcccwave.lasso.values <- matrix(NA, nrow = n, ncol = lasso.long)
for (j in 1:lasso.long) {
  lcccwave.lasso.values[, j] <- agree.cci(Y.mat, lcccwave.lasso[[1]][, , j])
}

wave.lasso.Chary_i_ <- apply(lcccwave.lasso.values, 2, mean)
wave.lasso.Chary1_i_ <- apply(lcccwave.lasso.values, 2, min)



texts1 <- expression(paste(bar(rho), " vs K"))
texts2 <- expression(paste(rho^0, " vs K"))


cex.lab <- 0.9 ### -> Size of labels for both x and y axis!
cex.axis <- 0.65 ### -> Size of coordinates for both x and y axis!
cex.main <- 0.7 ### -> Size of main topic ( USELESS since we will cut title )
cex.in <- 0.5 ### -> Size of all in boxplot (Ex: outlier )

yaxis.at <- c(-2, seq(0.85, 1.0, by = 0.05), 1.5)
yaxis.lab <- c("", "0.85", "0.90", "0.95", "1.00", "")
ylim <- c(0.84, 1.0)



itbasis <- sort(unique(c(lasso.x, lccc.pca[[2]])), decreasing = TRUE)

itxaxis <- seq(itbasis)
itidx1 <- rep(NA, length(itbasis))
itidx3 <- rep(NA, length(itbasis))
for (i in 1:length(itbasis)) {
  if (sum(lasso.x == itbasis[i]) != 0) {
    itidx1[i] <- i
  }
  if (sum(lccc.pca[[2]] == itbasis[i]) != 0) {
    itidx3[i] <- i
  }
}

xaxis.at <- c(-5, itxaxis[seq(3, length(itbasis), 3)], 35)
xaxis.lab <- c("", itbasis[ itxaxis[seq(3, length(itbasis), 3)] ], "")

xaxis.lab1 <- c("", lasso.x1[ itxaxis[seq(3, length(itbasis), 3)] ], "")


lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x1 <- lasso.counts.fit[[3]][ (lasso.counts.fit[[2]] < 1000) ]

Times.over.cv <- (Sys.time() - start)

########### start Reproduce Figure 3 ################################################################################


tiff("Figure3.tiff", width = 8, height = 4, units = "in", res = 120, bg = "transparent")

par(mfrow = c(1, 2), mar = c(4.5, 2.5, 3, 2))

plot(0, type = "o", xlab = "", ylab = "", cex.lab = 1.1, cex.axis = 0.7, axes = FALSE, ylim = ylim, xlim = c(0, length(itbasis) + 1))
points(itidx1[is.na(itidx1) != TRUE], lasso.Chary1_i_, col = "red", pch = 19)
lines(itidx1[is.na(itidx1) != TRUE], lasso.Chary1_i_, col = "red")
points(itidx1[is.na(itidx1) != TRUE], lasso.Chary_i_, col = "blue", pch = 19)
lines(itidx1[is.na(itidx1) != TRUE], lasso.Chary_i_, col = "blue")
axis(side = 2, at = yaxis.at, label = yaxis.lab, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
axis(side = 2, at = yaxis.at, label = rep("", length(yaxis.at)), line = -0.5, , tick = TRUE, las = 1)

axis(side = 1, at = xaxis.at, label = xaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis)

axis(side = 1, at = xaxis.at, label = rep("", length(xaxis.at)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at, label = xaxis.lab, line = 0.5, tick = FALSE, cex.axis = cex.axis)
mtext("C", side = 1, line = 0.5, at = 0.5, col = "blue", cex = 0.9)
mtext(expression(K[C]), side = 1, line = 1.7, at = 0.5, col = "blue", cex = 0.9)

title("A", cex.main = 1.2)


legend(4, 0.95, c(
  expression(paste(paste(paste(D^C, " ("), rho^0), ")")),
  expression(paste(paste(paste(D^C, " ("), bar(rho)), " )"))
),
lty = rep(1, 2), pch = rep(19, 2), col = c("red", "blue"),
cex = 1.1, bty = "n", ncol = 1
)


yaxis.at1 <- c(-2, seq(0.80, 1.0, by = 0.05), 1.5)
yaxis.lab1 <- c("", "0.80", "0.85", "0.90", "0.95", "1.00", "")
ylim1 <- c(0.80, 1.0)


plot(0, type = "o", xlab = "K", ylab = "", cex.lab = 1.1, cex.axis = 0.7, axes = FALSE, ylim = ylim1, xlim = c(0, length(itbasis) + 1))
points(seq(length(lccc.pca[[2]])), rev(apply(lccc.pca[[1]], 2, mean)), col = "gray", pch = 19)
lines(seq(length(lccc.pca[[2]])), rev(apply(lccc.pca[[1]], 2, mean)), col = "gray", lty = 2)
points(seq(length(lccc.pca[[2]])), rev(apply(lccc.pca[[1]], 2, min)), col = "gray", pch = 19)
lines(seq(length(lccc.pca[[2]])), rev(apply(lccc.pca[[1]], 2, min)), col = "gray")

points(itidx1, rev(wave.lasso.Chary1_i_), col = "red", pch = 19)
lines(itidx1, rev(wave.lasso.Chary1_i_), col = "red")
points(itidx1[is.na(itidx1) != TRUE], rev(wave.lasso.Chary_i_), col = "blue", pch = 19)
lines(itidx1[is.na(itidx1) != TRUE], rev(wave.lasso.Chary_i_), col = "blue")

axis(side = 2, at = yaxis.at1, label = yaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
axis(side = 2, at = yaxis.at1, label = rep("", length(yaxis.at1)), line = -0.5, , tick = TRUE, las = 1)

xaxis.at22 <- c(-5, itxaxis[ seq(1, length(itbasis), 3)], 35)
xaxis.lab22 <- c("", rev(itbasis)[  seq(1, length(itbasis), 3) ], "")

axis(side = 1, at = xaxis.at22, label = rep("", length(xaxis.at22)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at22, label = (xaxis.lab22), line = -0.5, tick = FALSE, cex.axis = cex.axis)

legend(6, 0.94, c(
  expression(paste("Quantlets (", rho^0, ")")),
  expression(paste("Quantlets (", bar(rho), " )")),
  expression(paste("PCA (", rho^0, ")")),
  expression(paste("PCA (", bar(rho), " )"))
),
lty = c(rep(1, 3), 2), pch = c(rep(19, 4)), col = c("red", "blue", "gray", "gray"),
cex = 1.1, bty = "n", ncol = 1
)
title("B", cex.main = 1.2)


dev.off()

########### end Reproduce Figure 4 ################################################################################



############# We choose # of basis ###############################################################################

lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x.idx <- seq(length(lasso.list1))[ (lasso.counts.fit[[2]] < 1000) ]
list.order <- lasso.counts.fit[[1]]
unlist(lapply(list.order, length))



be <- 3
REDUCED_BASE9 <- BETA_BASE_TOTAL_2[, list.order[[be + 8 ] ]] # our choice in paper
REDUCED_BASE21 <- BETA_BASE_TOTAL_2[, list.order[[be + 20 ] ]] # Normal case in paper


############ rearrange basis order #################################


Sys.time() -> start
grids.total.list <- vector("list", 64)
for (i in 1:n) {
  y <- raw.dataset[[i]]
  y.long <- length(y)
  grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
  grids.total.list[[i]] <- unique(round(grid.p, 4))
}

cbind(unlist(lapply(grids.total.list, length)), unlist(lapply(raw.dataset, length)), unlist(lapply(grids.total.list, length)) == unlist(lapply(raw.dataset, length)))
grids.total <- seq(0.0001, 0.9999, 0.0001)
grids.total.obs <- sort(unique(unlist(grids.total)))


### here we unify grid set as p1024, otherwise, it made too many problem. gramshumit, memory error, irregar, empiricalQ
grids.total <- p1024
CDFBETA <- GENERATE_BETA_CDF(a1, a2, grids.total)
NQ <- qnorm(grids.total)
BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))



############ compute quanvelts basis #################################

REDUCED_BASE9 <- BETA_BASE_TOTAL_2[, list.order[[be + 8 ] ]]

Gram9 <- gramSchmidt(REDUCED_BASE9, tol = .Machine$double.eps^0.5)
norms <- gramSchmidt(as.matrix(REDUCED_BASE21), tol = .Machine$double.eps^0.5)$Q
Quantlet9 <- WaveletDenose(Gram9$Q[, -1], filter.number = 2, family = "DaubExPhase")[, , 2]
Quantlet9s <- cbind(norms, centering.function(Quantlet9, scale = TRUE))
Quantlet21s <- norms


############ Reproduce Figure 4###############################################################################################

tiff("Figure4.tiff",
  width = 16, height = 12, units = "in", res = 160, bg = "transparent"
)
par(mfrow = c(4, 4), mar = c(4.5, 2.5, 3, 2))
plot(p1024, rep(1, 1024), type = "l", lty = 1, lwd = 0.2, main = bquote("Quantlet" ~ psi[.(1)]))
for (v in 1:15) {
  if (v >= 7) {
    ylims <- c(-0.05, 0.05)
  }
  if (v < 7) {
    ylims <- c(-0.2, 0.2)
  }
  plot(p1024, Quantlet9s[, v], type = "l", lty = 1, lwd = 0.2, main = bquote("Quantlet" ~ psi[.(v + 1)]), ylim = ylims, xlab = "")
}
dev.off()


################ compute empirical coefficients ###########################################################################


EmpCoefs_9 <- EmpCoefs(Quantlet9s, Y.list = raw.dataset1)

EmpCoefs <- function(B, Y.list) {
  n <- length(Y.list)
  mats <- matrix(NA, nrow = n, ncol = (dim(B)[2] + 1))
  for (i in 1:n) {
    y <- Y.list[[i]]
    y.long <- length(y)
    set.seed(123 + i)
    grids <- sort(sample(dim(B)[1], y.long))
    Psi <- cbind(rep(1, length(grids)), B[ grids, ])
    mats[i, ] <- solve(t(Psi) %*% Psi) %*% t(Psi) %*% y
  }
  return(mats)
}

Times.over.emp <- (Sys.time() - start)

########### import covariates and manage them to analysis ##########################################################################################


Sys.time() -> start
library(survival)

DATA <- read.csv("Q3_Covariates.csv", header = TRUE)

cut.time <- 12
stime <- as.numeric(DATA[, 6])
delta <- ifelse(as.character(DATA[, 7]) == "DECEASED", 1, 0)
kmfit <- survfit(Surv(stime, delta) ~ 1)
St <- summary(kmfit, time = c(cut.time, stime))$surv
Ft <- (1 - St)
up.delta <- ifelse((delta == 1) & (stime <= cut.time), 1, ifelse(stime > cut.time, 0, Ft[1] - Ft[-1]))
int <- rep(1, 64)
G <- ifelse(DATA$Gender == "FEMALE", 0, 1) ## female=0 male=1
A1 <- cbind(DATA$Classical, DATA$Mesenchymal, DATA$Neural, DATA$Proneural)
A2 <- cbind(DATA$DDIT3, DATA$EGFR, DATA$KIT, DATA$MDM4, DATA$PDGFRA, DATA$PIK3CA, DATA$PTEN)
AGE <- (DATA$Age - min(DATA$Age)) / (max(DATA$Age) - min(DATA$Age))
KFS <- (DATA$KF.score - min(DATA$KF.score)) / (max(DATA$KF.score) - min(DATA$KF.score))
LX <- cbind(int, G, A1, A2, (AGE - median(AGE)), (KFS - median(KFS)), (up.delta))
colnames(LX) <- c("i", "sex", "TC", "TM", "TN", "TP", "GD", "GE", "GK", "GM", "GPD", "GPI", "GPT", "age", "kpf", "del")
x <- LX[, c(1, 2, 4, 7, 8, 14, 16)]

X <- x

################## estimation #############################################################################################


Emp_fit_9 <- EmpQuant2(EmpCoefs_9, REDUCED_BASE9, Quantlet9s, X, delta2 = 0.95, H = 7)


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

Cluster_9 <- Cluster(REDUCED_BASE9, Gram9$Q, H = 7)



mcmc_fit_9 <- MCMC_QUANTLET(X, Y = Emp_fit_9$sd_l2, Cluster_9, Emp_fit_9$TB00, zeroKept = FALSE, n.iter = 2000, burn = 200)


mm_1 <- c(1, 1, 0.5, 0.5, 0.5, 0, 0.5)
mm_2 <- c(1, 0, 0.5, 0.5, 0.5, 0, 0.5)
mm_3 <- c(1, 0.5, 0.5, 0.5, 0.5, max(X[, 6]), 0.5)
mm_4 <- c(1, 0.5, 0.5, 0.5, 0.5, 0, 0.5)
mm_5 <- c(1, 0.5, 0.5, 1, 0.5, 0, 0.5)
mm_6 <- c(1, 0.5, 0.5, 0, 0.5, 0, 0.5)
mm_7 <- c(1, 0.5, 0.5, 0.5, 1, 0, 0.5)
mm_8 <- c(1, 0.5, 0.5, 0.5, 0, 0, 0.5)
mm_9 <- c(1, 0.5, 1, 0.5, 0.5, 0, 0.5)
mm_10 <- c(1, 0.5, 0, 0.5, 0.5, 0, 0.5)
mm_11 <- c(1, 0.5, 0.5, 0.5, 0.5, 0, 1)
mm_12 <- c(1, 0.5, 0.5, 0.5, 0.5, 0, 0)


PX0 <- rbind(mm_1, mm_2, mm_3, mm_4, mm_5, mm_6, mm_7, mm_8, mm_9, mm_10, mm_11, mm_12)
X1 <- PX0


mcmcInfer_9 <- inference(mcmc_fit_9[[1]], BackTransfor = Emp_fit_9$sdPhi, X, signifit = 0.975, X1 = PX0, p = p1024, n.sup = 77, xranges = c(-10, 60))


Ip <- diag(rep(1, length(p1024)), length(p1024), length(p1024))


##  MCMC_PCR  -- conducts MCMC computaions based on PC regression

mcmc_fit_22 <- MCMC_PCR(X, t(Qy), BackTransfor = Ip, n.iter = 2000, burn = 200)
mcmcInfer_22 <- inference(mcmc_fit_22[[1]], BackTransfor = Ip, X, signifit = 0.975, X1 = PX0, p = p1024, n.sup = 77, xranges = c(-10, 60))


mcmc_fit_9_norm <- MCMC_NRPCT(mcmc_fit_9[[1]], 0.975, PX0)

mcmc_fit_22_norm <- MCMC_NRPCT(mcmc_fit_22[[1]], 0.975, PX0)


Times.over.est <- (Sys.time() - start)



source("plots.R")

n.sup <- 77
xdomain <- seq(-10, 60, length.out = n.sup)

############ Reproduce Figure 6 ###############################################################################################
tiff(
  "Figure6.tiff",
  width = 16, height = 12, units = "in", res = 150, bg = "transparent"
)
RealPlots(mcmcInfer_9, p1024, edit = 15, opt = 1)
dev.off()

############ Reproduce Figure 7 ###############################################################################################
tiff(
  "Figure7.tiff",
  width = 8, height = 5, units = "in", res = 150, bg = "transparent"
)

opt <- 1

par(mfrow = c(1, 2), mar = c(2.5, 2.5, 3, 2))

plot(0, type = "n", ylim = c(-15, 45), xlim = c(0, 1), main = "")
lines(p1024, mcmcInfer_9$DataEst[, 4], col = "hotpink", lty = 1, lwd = 2)
lines(p1024, mcmcInfer_9$estCIu[, 4], col = "gray", lty = 2, lwd = 2)
lines(p1024, mcmcInfer_9$estCIl[, 4], col = "gray", lty = 2, lwd = 2)
lines(p1024, mcmcInfer_9$jointCI[, 4, 1 ], col = "black", lty = 3, lwd = 2)
lines(p1024, mcmcInfer_9$jointCI[, 4, 2 ], col = "black", lty = 3, lwd = 2)

title("A", cex = 1.5)
points(p1024[mcmcInfer_9$local_p[, 3]], -15 * mcmcInfer_9$local_p[, 3][mcmcInfer_9$local_p[, 3] == TRUE], col = "orange")


legend(0.1, 47,
  pch = c(NA, NA, NA, 19),
  c("DDIT3 effect", "95% CI (point)", "95% CI (joint)", "Flag"),
  lty = c(1, 2, 3, 1),
  col = c("hotpink", "gray", "black", "orange"),
  cex = 1.5, bty = "n", ncol = 1
)

plot(0, type = "n", ylim = c(-15, 45), xlim = c(0, 1), main = "")
lines(p1024, mcmcInfer_22$DataEst[, 4], col = "hotpink", lty = 1, lwd = 2)
lines(p1024, mcmcInfer_22$estCIu[, 4], col = "gray", lty = 2, lwd = 2)
lines(p1024, mcmcInfer_22$estCIl[, 4], col = "gray", lty = 2, lwd = 2)
lines(p1024, mcmcInfer_22$jointCI[, 4, 1 ], col = "black", lty = 3, lwd = 2)
lines(p1024, mcmcInfer_22$jointCI[, 4, 2 ], col = "black", lty = 3, lwd = 2)

title("B", cex = 1.5)
points(p1024[mcmcInfer_22$local_p[, 3]], -15 * mcmcInfer_22$local_p[, 3][mcmcInfer_22$local_p[, 3] == TRUE], col = "orange")

legend(0.1, 47,
  pch = c(NA, NA, NA, NA),
  c("DDIT3 effect", "95% CI (point)", "95% CI (joint)", NA),
  lty = c(1, 2, 3, NA),
  col = c("hotpink", "gray", "black", NA),
  cex = 1.5, bty = "n", ncol = 1
)

dev.off()
################### conduct feather selection method ##########################################





################### conduct feather selection method ##########################################
library(moments)

quantiles2 <- function(x, probs = c(0.025, 0.975)) {
  quantile(x, probs)
}

feature_mean <- apply(Qy, 2, mean)
feature_var <- apply(Qy, 2, var)^{
  0.5
}
feature_skew <- apply(Qy, 2, skewness)
feature_kurto <- apply(Qy, 2, kurtosis)
feature_stats <- apply(Qy, 2, summary)[ c(2, 3, 5), ]

Features <- cbind(feature_mean, feature_var, feature_skew, feature_kurto, t(feature_stats))

n.iter <- 2000
N <- 64
ITF_UI_FEATURE_T <- array(NA, c(n.iter, dim(Features)[2], dim(X)[2]))

nu0 <- 0.0064
s20 <- 1500
Px <- dim(X)[2]

for (j in 1:dim(Features)[2]) { ##  j =1
  set.seed(3181 + j)
  f_j <- as.numeric(Features[, j], ncol = 1)
  Hg <- X %*% solve(t(X) %*% X) %*% t(X)
  SSRgf <- t(f_j) %*% (diag(1, nrow = N) - Hg) %*% f_j
  s2_f <- 1 / rgamma(n.iter, (nu0 + N) / 2, (0.001 + SSRgf) / 2)
  Vb <- solve(t(X) %*% X)
  Eb_f <- Vb %*% t(X) %*% f_j
  E_f <- matrix(rnorm(n.iter * Px, 0, sqrt(s2_f)), n.iter, Px)
  beta_f <- t(t(E_f %*% chol(Vb)) + c(Eb_f)) ## dim(beta)
  for (k in 1:Px) {
    ITF_UI_FEATURE_T[, j, k] <- beta_f[, k]
  }
}

burn <- 200
ITF_UI_FEATURE <- ITF_UI_FEATURE_T[-c(1:burn), , ]

FEATURE_MU <- ITF_UI_FEATURE[, 1, ] %*% t(PX0)
STAT_MU <- array(NA, c((n.iter - burn), dim(PX0)[1], 1))
for (i in 1:(n.iter - burn)) {
  STAT_MU[i, , 1 ] <- FEATURE_MU[i, ]
}
MU_FTHR <- rbind(apply(STAT_MU[, , 1], 2, mean), apply(STAT_MU[, , 1], 2, quantiles2))

FEATURE_VAR <- ITF_UI_FEATURE[, 2, ] %*% t(PX0)
STAT_VAR <- array(NA, c((n.iter - burn), dim(PX0)[1], 1))
for (i in 1:(n.iter - burn)) {
  STAT_VAR[i, , 1  ] <- FEATURE_VAR[i, ]
}
VAR_FTHR <- rbind(apply(STAT_VAR[, , 1], 2, mean), apply(STAT_VAR[, , 1], 2, quantiles2))

FEATURE_MU3 <- ITF_UI_FEATURE[, 3, ] %*% t(PX0)
STAT_MU3 <- array(NA, c((n.iter - burn), dim(PX0)[1], 4))
for (i in 1:(n.iter - burn)) {
  STAT_MU3[i, , 1  ] <- FEATURE_MU3[i, ]
}
MU3_FTHR <- rbind(apply(STAT_MU3[, , 1], 2, mean), apply(STAT_MU3[, , 1], 2, quantiles2))


spx <- dim(PX0)[1]
STAT_MU_DIFF <- matrix(NA, nrow = spx / 2, ncol = 9)

TRUE_MU <- g.mean
k.n <- 0
for (k in 1:(spx / 2)) {
  i <- 2 * (k - 1) + 1
  j <- 2 * (k - 1) + 2
  k.n <- k.n + 1
  TEST_STATS <- STAT_MU[, j, 1 ] - STAT_MU[, i, 1]


  M_j <- mean(TEST_STATS)
  S_j <- var(TEST_STATS)^{
    0.5
  }
  Z_j <- (TEST_STATS - M_j) / S_j
  Quan_Z_j_u <- quantile(Z_j, 0.975)
  Quan_Z_j_l <- quantile(Z_j, 0.025)

  STAT_MU_DIFF[k.n, 1] <- i
  STAT_MU_DIFF[k.n, 2] <- j
  STAT_MU_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
  STAT_MU_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
  STAT_MU_DIFF[k.n, 8] <- 2 * mean(abs(M_j / S_j) <= Z_j)
  STAT_MU_DIFF[k.n, 9] <- ifelse(STAT_MU_DIFF[k.n, 3] == 0, 0, 1)
}

STAT_VAR_DIFF <- matrix(NA, nrow = spx / 2, ncol = 9)
k.n <- 0
for (k in 1:(spx / 2)) {
  i <- 2 * (k - 1) + 1
  j <- 2 * (k - 1) + 2
  k.n <- k.n + 1
  TEST_STATS <- STAT_VAR[, j, 1] - STAT_VAR[, i, 1]

  M_j <- mean(TEST_STATS)
  S_j <- var(TEST_STATS)^{
    0.5
  }
  Z_j <- (TEST_STATS - M_j) / S_j
  Quan_Z_j_u <- quantile(Z_j, 0.975)
  Quan_Z_j_l <- quantile(Z_j, 0.025)

  STAT_VAR_DIFF[k.n, 1] <- i
  STAT_VAR_DIFF[k.n, 2] <- j
  STAT_VAR_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
  STAT_VAR_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
  STAT_VAR_DIFF[k.n, 8] <- 2 * mean(abs(M_j / S_j) <= Z_j)
  STAT_VAR_DIFF[k.n, 9] <- ifelse(STAT_VAR_DIFF[k.n, 3] == 0, 0, 1)
}


STAT_MU3_DIFF <- matrix(NA, nrow = spx / 2, ncol = 9)
k.n <- 0
for (k in 1:(spx / 2)) {
  i <- 2 * (k - 1) + 1
  j <- 2 * (k - 1) + 2
  k.n <- k.n + 1
  TEST_STATS <- STAT_MU3[, j, 1] - STAT_MU3[, i, 1 ]

  M_j <- mean(TEST_STATS)
  S_j <- var(TEST_STATS)^{
    0.5
  }
  Z_j <- (TEST_STATS - M_j) / S_j
  Quan_Z_j_u <- quantile(Z_j, 0.975)
  Quan_Z_j_l <- quantile(Z_j, 0.025)

  STAT_MU3_DIFF[k.n, 1] <- i
  STAT_MU3_DIFF[k.n, 2] <- j
  STAT_MU3_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
  STAT_MU3_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
  STAT_MU3_DIFF[k.n, 8] <- 2 * mean(abs(M_j / S_j) <= Z_j)
  STAT_MU3_DIFF[k.n, 9] <- ifelse(STAT_MU3_DIFF[k.n, 3] == 0, 0, 1)
}


################### reproduce table 4 ##########################################


table4 <- round(
  cbind(
    cbind(mcmcInfer_22$mu_diff[, 8], mcmcInfer_9$mu_diff[, 8], STAT_MU_DIFF[, 8]),

    cbind(mcmcInfer_22$sigma_diff[, 8], mcmcInfer_9$sigma_diff[, 8], STAT_VAR_DIFF[, 8]),

    cbind(mcmcInfer_22$mu3_diff[, 8], mcmcInfer_9$mu3_diff[, 8], STAT_MU3_DIFF[, 8])
  ), 3
)

colnames(table4) <- c("B_mu", "E_mu", "G_mu", "B_se", "E_se", "G_se", "B_mu3", "E_mu3", "G_mu3")
rownames(table4) <- c("Sex", "Age", "DDIT", "EGFR", "Mesen", "Surv")

write.csv(table4, "Table4.csv")
