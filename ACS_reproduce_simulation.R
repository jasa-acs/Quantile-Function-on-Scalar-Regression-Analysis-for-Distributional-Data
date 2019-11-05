
rm(list = ls())

library(sn)
library(glmnet)
library(mvtnorm)
library(pracma)
library(matrixcalc)
library(MASS)

library(pracma)
library(wavethresh)



source("getquantlets.R")
source("qfreg.R")
source("inference.R")
source("PcrQuant.R")
source("datagen.R")

###################################################################################################

p1024 <- signif(seq(0.001, 0.999, length = 1024), 4)
p <- p1024
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

BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))

BETASCDF1 <- t(GENERATE_BETA_CDF(a1, a2, p1024))
BETA_BASE_TOTAL_2 <- cbind(BNQ, BETASCDF1)


raw.dataset <- vector("list", n)
for (i in 1:n) {
  raw.dataset[[i]] <- QY[, i]
}



########################  lasso selection! ########################

n <- 120
Sys.time() -> start
lasso.list <- vector("list", n)
for (i in 1:n) {
  y <- raw.dataset[[i]]
  y.long <- length(y)
  set.seed(12345 + i)
  lasso_fit <- glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE)
  cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3)
  zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
  selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros == FALSE ]
  lasso.list[[i]] <- selects
}
Times.over.lasso <- (Sys.time() - start)




########### ################################################################################################

Sys.time() -> start
lasso.list1 <- lapply(lasso.list, CatNorm)
lasso.nonzero.obs1 <- lapply(lasso.list1, rep.list)

lasso.counts.fit <- CountBasis(lasso.list1, lasso.nonzero.obs1)

n <- 120


lasso_IncidenceVec_i_ <- vector("list", n)
for (i in 1:n) {
  lasso_fitIncd <- IncidenceVec(lasso.list1[-i], lasso.nonzero.obs1[-i])
  lasso_IncidenceVec_i_[[i]] <- lasso_fitIncd
}



lccc.lasso <- LCCC.FIT(lasso_IncidenceVec_i_,
  B = BETA_BASE_TOTAL_2, remain.counts = lasso.counts.fit[[3]],
  remain.basis = lasso.counts.fit[[2]], raw.dataset
)

Qy <- QY

fit.pca <- TruncPCA(Qy, length.train = NULL, max.K = n)

lccc.pca <- LCCC.PCR.FIT(Y = Qy, B = fit.pca$v, remain.basis = lasso.counts.fit[[2]])



Y.mat <- matrix(NA, nrow = 1024, ncol = n)
for (i in 1:n) {
  Y.mat[(1:length(raw.dataset[[i]])), i] <- raw.dataset[[i]]
}


lasso.long <- sum(lasso.counts.fit[[2]] < 1000)

lasso.values <- matrix(NA, nrow = n, ncol = lasso.long)

for (j in 1:lasso.long) {
  lasso.values[, j] <- agree.cci(Y.mat, lccc.lasso[[1]][, , j])
}

lasso.Chary_i_ <- apply(lasso.values, 2, mean)
lasso.Chary1_i_ <- apply(lasso.values, 2, min)

lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x1 <- lasso.counts.fit[[3]][ (lasso.counts.fit[[2]] < 1000) ]



PROJECT <- "S2_"
texts1 <- expression(paste(bar(rho), " vs K"))
texts2 <- expression(paste(rho^0, " vs K"))


cex.lab <- 0.9 ### -> Size of labels for both x and y axis!
cex.axis <- 0.65 ### -> Size of coordinates for both x and y axis!
cex.main <- 0.7 ### -> Size of main topic ( USELESS since we will cut title )
cex.in <- 0.5 ### -> Size of all in boxplot (Ex: outlier )

yaxis.at <- c(-2, seq(0.97, 1.0, by = 0.01), 1.5)
yaxis.lab <- c("", "0.97", "0.98", "0.99", "1.00", "")
ylim <- c(0.968, 1.0)



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




par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 2))

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

legend(4, 0.99, c(expression(rho^0), expression(bar(rho))),
  lty = rep(1, 2), pch = rep(19, 2), col = c("red", "blue"),
  cex = 1.1, bty = "n", ncol = 1
)





###################################################################################################





lcccwave.lasso <- LCCC.WAVE.FIT(
  global.list = lasso.counts.fit[[1]], B = BETA_BASE_TOTAL_2,
  remain.counts = lasso.counts.fit[[3]],
  remain.basis = lasso.counts.fit[[2]], Y.list = raw.dataset
)


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

yaxis.at1 <- c(-2, seq(0.97, 1.0, by = 0.01), 1.5)
yaxis.lab1 <- c("", "0.97", "0.98", "0.99", "1.00", "")
ylim1 <- c(0.968, 1.0)



par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 2))

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

xaxis.at22 <- c(-5, itxaxis[ seq(2, length(itbasis), 3)], 35)
xaxis.lab22 <- c("", rev(itbasis)[  seq(2, length(itbasis), 3) ], "")

axis(side = 1, at = xaxis.at22, label = rep("", length(xaxis.at22)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at22, label = (xaxis.lab22), line = -0.5, tick = FALSE, cex.axis = cex.axis)

legend(6, 0.99, c(
  expression(paste("Quantlets (", rho^0, ")")),
  expression(paste("Quantlets (", bar(rho), " )")),
  expression(paste("PCA (", rho^0, ")")),
  expression(paste("PCA (", bar(rho), " )"))
),
lty = c(rep(1, 3), 2), pch = c(rep(19, 4)), col = c("red", "blue", "gray", "gray"),
cex = 1.1, bty = "n", ncol = 1
)





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


legend(4, 0.99, c(
  expression(paste(paste(paste(D^u, " ("), rho^0), ")")),
  expression(paste(paste(paste(D^u, " ("), bar(rho)), " )"))
),
lty = rep(1, 2), pch = rep(19, 2), col = c("red", "blue"),
cex = 1.1, bty = "n", ncol = 1
)



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

xaxis.at22 <- c(-5, itxaxis[ seq(2, length(itbasis), 3)], 35)
xaxis.lab22 <- c("", rev(itbasis)[  seq(2, length(itbasis), 3) ], "")

axis(side = 1, at = xaxis.at22, label = rep("", length(xaxis.at22)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at22, label = (xaxis.lab22), line = -0.5, tick = FALSE, cex.axis = cex.axis)

legend(6, 0.99, c(
  expression(paste("Quantlets (", rho^0, ")")),
  expression(paste("Quantlets (", bar(rho), " )")),
  expression(paste("PCA (", rho^0, ")")),
  expression(paste("PCA (", bar(rho), " )"))
),
lty = c(rep(1, 3), 2), pch = c(rep(19, 4)), col = c("red", "blue", "gray", "gray"),
cex = 1.1, bty = "n", ncol = 1
)
title("B", cex.main = 1.2)





Times.over.cv <- (Sys.time() - start)

#########################################################################################################


Sys.time() -> start
lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x.idx <- seq(length(lasso.list1))[ (lasso.counts.fit[[2]] < 1000) ]
list.order <- lasso.counts.fit[[1]]
unlist(lapply(list.order, length))


be <- 2


REDUCED_BASE18 <- BETA_BASE_TOTAL_2[, list.order[[be + 17 ] ]]
REDUCED_BASE19 <- BETA_BASE_TOTAL_2[, list.order[[be + 18 ] ]]

REDUCED_BASE28 <- BETA_BASE_TOTAL_2[, list.order[[be + 27 ] ]]




list.order.up <- vector("list", 28)
for (i in 3:25) { 
  emp0 <- EmpCoefs(BETA_BASE_TOTAL_2[, list.order[[ be + i ] ] ], Y.list = raw.dataset)
  list.order.up[[(i + be)]] <- uporder(emp0)
}
list.order.up[[28]] <- list.order[[ be + 27 ] ] ## NA is becuase from singular


REDUCED_BASE18 <- BETA_BASE_TOTAL_2[, list.order.up[[18 ] ]]
REDUCED_BASE19 <- BETA_BASE_TOTAL_2[, list.order.up[[19 ] ]]


Gram18 <- gramSchmidt(REDUCED_BASE18, tol = .Machine$double.eps^0.5)
Gram28 <- gramSchmidt(as.matrix(REDUCED_BASE28), tol = .Machine$double.eps^0.5)


Quantlet18 <- WaveletDenose(Gram18$Q[, -1], filter.number = 2, family = "DaubExPhase")[, , 2]
norms <- Gram28$Q


Quantlet18s <- cbind(norms, centering.function(Quantlet18, scale = TRUE))
Quantlet28s <- norms

Times.over.emp <- (Sys.time() - start)






Sys.time() -> start
EmpCoefs_18 <- EmpCoefs(Quantlet18s, Y.list = raw.dataset)
EmpCoefs_28 <- EmpCoefs(Quantlet28s, Y.list = raw.dataset)


X <- x
X1 <- x[c(1, 61, 31, 91), ]
PX0 <- X1



Emp_fit_18 <- EmpQuant2(EmpCoefs_18, REDUCED_BASE18, Quantlet18s, X, delta2 = 0.95, H = 7)
Emp_fit_28 <- EmpQuant2(EmpCoefs_28, as.matrix(REDUCED_BASE28), as.matrix(Quantlet28s), X, delta2 = 1.95, H = 7)


mcmc_fit_18 <- MCMC_QUANTLET(X, Y = Emp_fit_18$sd_l2, Emp_fit_18$cluster, Emp_fit_18$TB00, zeroKept = FALSE, n.iter = 2000, burn = 200)
mcmc_fit_28 <- MCMC_PCR(X, Emp_fit_28$sd_l2, Emp_fit_28$sdPhi, n.iter = 2000, burn = 200)



mcmcInfer_18 <- inference(mcmc_fit_18[[1]], BackTransfor = Emp_fit_18$sdPhi, X, signifit = 0.975, X1 = PX0, p = p1024, n.sup = 53, xranges = c(-30, 30))
mcmcInfer_28 <- inference(mcmc_fit_28[[1]], BackTransfor = Emp_fit_28$sdPhi, X, signifit = 0.975, X1 = PX0, p = p1024, n.sup = 53, xranges = c(-30, 30))



#########################################################################################################





Truns <- c(20)
Phi1 <- fit.pca$v[, 1:Truns[1] ]
scale1 <- diag(1 / fit.pca$d[1:Truns[1]]^0.5, Truns[1], Truns[1]) / sqrt(dim(Qy)[2] - 1)
sPhi1 <- Phi1 %*% scale1

Qy_mean <- fit.pca$mu
Qy_mean_mat <- matrix(rep(Qy_mean, each = dim(Qy)[2]), nrow = dim(Qy)[1], ncol = dim(Qy)[2], byrow = TRUE)
Res_Qy <- as.matrix(Qy) - Qy_mean_mat

Z1 <- t(Res_Qy) %*% Phi1 
U1 <- t(Res_Qy) %*% sPhi1 
int <- rep(1, dim(Qy)[2])
Zi1 <- cbind(int, Z1) 
Ui1 <- cbind(int, U1)
Ip <- diag(rep(1, length(p)), length(p), length(p))


mcmc_fit_29 <- MCMC_PCR(X, t(Qy), BackTransfor = Ip, n.iter = 2000, burn = 200)
mcmcInfer_29 <- inference(mcmc_fit_29[[1]], BackTransfor = Ip, X, signifit = 0.975, X1 = PX0, p = p1024, n.sup = 53, xranges = c(-30, 30))


mcmc_fit_30 <- MCMC_PCR(X, Ui1, cbind(Qy_mean, Phi1 %*% solve(scale1)), n.iter = 2000, burn = 200)
mcmcInfer_30 <- inference(mcmc_fit_30[[1]], BackTransfor = cbind(Qy_mean, Phi1 %*% solve(scale1)), X, signifit = 0.975, X1 = PX0, p = p1024, n.sup = 53, xranges = c(-30, 30))


mcmc_fit_32 <- MCMC_PCR(X, Emp_fit_18$sd_l2, Emp_fit_18$sdPhi, n.iter = 2000, burn = 200)
mcmcInfer_32 <- inference(mcmc_fit_32[[1]], BackTransfor = Emp_fit_18$sdPhi, X, signifit = 0.975, X1 = PX0, p = p1024, n.sup = 53, xranges = c(-30, 30))



quantMCMC <- MCMC_QR(X = x, y = QY, n.iter = 500, burn = 200, ep1 = 0.0064, tau = 0.5)
mcmcInfer_34 <- inference(quantMCMC[[1]], BackTransfor = Ip, X, signifit = 0.975, X1 = PX0, p = p1024, n.sup = 53, xranges = c(-30, 30))


Times.over.est <- (Sys.time() - start)
##############################################################################################################



b1 <- qdf1
b2 <- qdf2 - qdf1
b3 <- qdf3 - qdf1
b4 <- qdf4 - qdf1

ss <- c(34, 29, 30, 32, 18, 28)


Cvrg_fit <- matrix(NA, ncol = 6, nrow = 4)
Area_fit <- matrix(NA, ncol = 6, nrow = 4)


for (i in 1:length(ss)) { 
  at <- ss[i]
  obj_name <- paste0("mcmcInfer_", at, sep = "")
  obj_output <- eval(parse(text = obj_name))
  Cvrg_fit[, i] <- c(
    sum(ifelse(b1 >= obj_output$jointCI[, 1, 2], 1, 0) * ifelse(b1 <= obj_output$jointCI[, 1, 1], 1, 0)) / length(b1),
    sum(ifelse(b2 >= obj_output$jointCI[, 2, 2], 1, 0) * ifelse(b2 <= obj_output$jointCI[, 2, 1], 1, 0)) / length(b1),
    sum(ifelse(b3 >= obj_output$jointCI[, 3, 2], 1, 0) * ifelse(b3 <= obj_output$jointCI[, 3, 1], 1, 0)) / length(b1),
    sum(ifelse(b4 >= obj_output$jointCI[, 4, 2], 1, 0) * ifelse(b4 <= obj_output$jointCI[, 4, 1], 1, 0)) / length(b1)
  )
  Area_fit[, i] <- apply(obj_output$jointCI[, , 1] - obj_output$jointCI[, , 2], 2, mean)
} 


table2_1 <- round(Area_fit, 3)
colnames(table2_1) <- c("A", "B", "C", "D", "E", "F")
rownames(table2_1) <- c("beta_1", "beta_2", "beta_3", "beta_4")
write.csv(table2_1, "Table2_Area.csv")

table2_2 <- round(Cvrg_fit, 3)
colnames(table2_2) <- c("A", "B", "C", "D", "E", "F")
rownames(table2_2) <- c("beta_1", "beta_2", "beta_3", "beta_4")
write.csv(table2_2, "Table2_coverage.csv")






library(moments)

feature_mean <- apply(Qy, 2, mean)
feature_var <- apply(Qy, 2, var)^{
  0.5
}
feature_skew <- apply(Qy, 2, skewness)
feature_kurto <- apply(Qy, 2, kurtosis)
feature_stats <- apply(Qy, 2, summary)[ c(2, 3, 5), ]

Features <- cbind(feature_mean, feature_var, feature_skew, feature_kurto, t(feature_stats))

###########################################################################################################
n.iter <- 2000
N <- 120
ITF_UI_FEATURE_T <- array(NA, c(n.iter, dim(Features)[2], dim(X)[2])) 

nu0 <- 0.0064
s20 <- 1500
Px <- 4

for (j in 1:dim(Features)[2]) { 
  set.seed(3181 + j)
  f_j <- as.numeric(Features[, j], ncol = 1)
  Hg <- X %*% solve(t(X) %*% X) %*% t(X)
  SSRgf <- t(f_j) %*% (diag(1, nrow = N) - Hg) %*% f_j
  s2_f <- 1 / rgamma(n.iter, (nu0 + N) / 2, (0.001 + SSRgf) / 2)
  Vb <- solve(t(X) %*% X)
  Eb_f <- Vb %*% t(X) %*% f_j
  E_f <- matrix(rnorm(n.iter * Px, 0, sqrt(s2_f)), n.iter, Px)
  beta_f <- t(t(E_f %*% chol(Vb)) + c(Eb_f)) 
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
  STAT_MU_DIFF[k.n, 3] <- round(TRUE_MU[j] - TRUE_MU[i], 2)
  STAT_MU_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
  STAT_MU_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
  STAT_MU_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
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
  ## STAT_VAR_DIFF[k.n, 3]= round( TRUE_VAR[j] -  TRUE_VAR[i] , 2)
  STAT_VAR_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
  STAT_VAR_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
  STAT_VAR_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
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
  ## STAT_MU3_DIFF[k.n, 3]= round( TRUE_MU3[j] -  TRUE_MU3[i] , 2)
  STAT_MU3_DIFF[k.n, 4:6] <- c(M_j, M_j + Quan_Z_j_l * S_j, M_j + Quan_Z_j_u * S_j)
  STAT_MU3_DIFF[k.n, 7] <- ifelse((M_j + Quan_Z_j_l * S_j < 0) & (M_j + Quan_Z_j_u * S_j > 0) == TRUE, 0, 1)
  STAT_MU3_DIFF[k.n, 8] <- mean(abs(M_j / S_j) <= Z_j)
  STAT_MU3_DIFF[k.n, 9] <- ifelse(STAT_MU3_DIFF[k.n, 3] == 0, 0, 1)
}

################### reproduce table 4 ##########################################


ss <- c(34, 29, 30, 32, 18, 28)




table3 <- round(
  cbind(
    c(mcmcInfer_34$mu_diff[, 8], mcmcInfer_34$sigma_diff[, 8], mcmcInfer_34$mu3_diff[, 8]),

    c(mcmcInfer_29$mu_diff[, 8], mcmcInfer_29$sigma_diff[, 8], mcmcInfer_29$mu3_diff[, 8]),
    c(mcmcInfer_30$mu_diff[, 8], mcmcInfer_30$sigma_diff[, 8], mcmcInfer_30$mu3_diff[, 8]),

    c(mcmcInfer_32$mu_diff[, 8], mcmcInfer_32$sigma_diff[, 8], mcmcInfer_32$mu3_diff[, 8]),

    c(mcmcInfer_18$mu_diff[, 8], mcmcInfer_18$sigma_diff[, 8], mcmcInfer_18$mu3_diff[, 8]),
    c(mcmcInfer_28$mu_diff[, 8], mcmcInfer_28$sigma_diff[, 8], mcmcInfer_28$mu3_diff[, 8]),
    c(STAT_MU_DIFF[, 8], STAT_VAR_DIFF[, 8], STAT_MU3_DIFF[, 8])
  ), 3
)


colnames(table3) <- c("A", "B", "C", "D", "E", "F", "G")
rownames(table3) <- c("m1=m3", "m2=m4", "s1=s3", "s2=s4", "x1=x3", "x2=x4")

write.csv(table3, "Table3.csv")




###############################################################################################################################





tiff(
  "Figure5.tiff",
  width = 12, height = 8, units = "in", res = 120, bg = "transparent"
)



par(mfrow = c(2, 3), mar = c(2.5, 2.5, 3, 2))

plot(0, xlim = c(min(z), 25), ylim = c(0, 0.15), type = "n", main = "A")
lines(z, pdf1, lty = 1, lwd = 2)
lines(z, pdf2, lty = 1, col = "red")
lines(z, pdf3, lty = 1, col = "blue")
lines(z, pdf4, lty = 1, col = "green")

legend(-20, 0.15, c("Normal", "Loc-shift", "Scale-diff", "Skew to L"),
  lty = rep(1, 4),
  col = c("black", "red", "blue", "green"),
  cex = 1.5, bty = "n", ncol = 1
)

title("A", cex.main = 1.2)


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

title("B", cex.main = 1.2)


legend(4, 0.99, c(
  expression(paste(paste(paste(D^C, " ("), rho^0), ")")),
  expression(paste(paste(paste(D^C, " ("), bar(rho)), " )"))
),
lty = rep(1, 2), pch = rep(19, 2), col = c("red", "blue"),
cex = 1.5, bty = "n", ncol = 1
)




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

xaxis.at22 <- c(-5, itxaxis[ seq(2, length(itbasis), 3)], 35)
xaxis.lab22 <- c("", rev(itbasis)[  seq(2, length(itbasis), 3) ], "")

axis(side = 1, at = xaxis.at22, label = rep("", length(xaxis.at22)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at22, label = (xaxis.lab22), line = -0.5, tick = FALSE, cex.axis = cex.axis)

legend(6, 0.99, c(
  expression(paste("Quantlets (", rho^0, ")")),
  expression(paste("Quantlets (", bar(rho), " )")),
  expression(paste("PCA (", rho^0, ")")),
  expression(paste("PCA (", bar(rho), " )"))
),
lty = c(rep(1, 3), 2), pch = c(rep(19, 4)), col = c("red", "blue", "gray", "gray"),
cex = 1.5, bty = "n", ncol = 1
)
title("C", cex.main = 1.2)


Qfit <- Emp_fit_18$sdPhi %*% t(Emp_fit_18$sd_l2)
ran <- sample(n, 45)
plot(0, type = "n", col = "orange", ylim = c(-15, 17), xlim = c(0, 1), lty = 3, main = "D")
for (i in 1:length(ran)) {
  lines(p, QY[, ran[i] ], col = "orange", lty = 1)
} 
for (i in 1:length(ran)) {
  lines(p, Qfit[, ran[i]], col = "black", lty = 3)
}


legend(0, 17, c("Obs Quantiles", "Quantlets"),
  lty = c(1, 3),
  col = c("orange", "black"),
  cex = 1.5, bty = "n", ncol = 1
)





plot(0, type = "n", col = "orange", ylim = c(-15, 15), xlim = c(0, 1), lty = 3, main = "E")
for (i in 1:n) {
  lines(p, QY[, i], col = "gray", lty = 3)
} 
lines(p, qdf1, col = "black", lty = 1, lwd = 2)
lines(p, qdf2, col = "red", lty = 1, lwd = 2)
lines(p, qdf3, col = "blue", lty = 1, lwd = 2)
lines(p, qdf4, col = "lightgreen", lty = 1, lwd = 2)

legend(0, 17, c("Overall quantile", "Loc-shift quantile", "Scale effect quantile", "Skew quantile"),
  lty = rep(1, 4),
  col = c("black", "red", "blue", "lightgreen"),
  cex = 1.5, bty = "n", ncol = 1
)

beta1 <- expression(paste(beta[1], ":overall effect"))
beta2 <- expression(paste(beta[2], ":location effect"))
beta3 <- expression(paste(beta[3], ":scale effect"))
beta4 <- expression(paste(beta[4], ":skew effect"))


plot(0, type = "n", col = "orange", ylim = c(-8, 10), xlim = c(0, 1), lty = 3, main = "F")
for (i in 1:n) {
  lines(p, QY[, i], col = "gray", lty = 3)
} 
lines(p, qdf1, col = "black", lty = 1, lwd = 2)
lines(p, qdf2 - qdf1, col = "red", lty = 1, lwd = 2)
lines(p, qdf3 - qdf1, col = "blue", lty = 1, lwd = 2)
lines(p, qdf4 - qdf1, col = "lightgreen", lty = 1, lwd = 2)

legend(-0.05, 11, c(beta1, beta2, beta3, beta4),
  lty = rep(1, 4),
  col = c("black", "red", "blue", "lightgreen"),
  cex = 1.5, bty = "n", ncol = 1
)



dev.off()


#####################################################################################
