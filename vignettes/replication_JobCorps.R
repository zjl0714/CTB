######################################################
############# Replication of JobCorps ################
######################################################

rm(list=ls())
library(dplyr)
library(grf)
library(hdm)
library(ggplot2)
library(sandwich)

setwd("/Users/yewang/Documents/GitHub/trimming/replication/")
source("application/CTTB.R")

seed <- 2023


load("application/data/JobCorps_baseline.rda")
load("application/data/JobCorps_employment.rda")
load("application/data/JobCorps_wages.rda")

# week 90
leedata=data.frame(treat=JobCorps_data_baseline$TREATMNT.y,
                   # selection=JobCorps_data_employment$week_208,
                   # outcome=JobCorps_data_wages$week_208,
                   selection=JobCorps_data_employment$week_90,
                   outcome=JobCorps_data_wages$week_90,
                   JobCorps_data_baseline[, 3:29])
leedata$logwage <- log(leedata$WKEARNR+1)

Y <- "outcome"
D <- "treat"
S <- "selection"
X <- names(leedata)[c(4:29, 31)]
Pscore <- "Pscore"

dat <- leedata[complete.cases(leedata[, D]), ]
dat <- dat[dat[, Y] != Inf, ]
for (i in 1:length(X)){
  Xi <- X[i]
  dat[, Xi] <- as.numeric(dat[, Xi])
}

Ntr <- sum(dat[, D] == 1)
N <- nrow(dat)

dat$Pscore <- Ntr/N

pf0 <- probability_forest(X = as.matrix(dat[dat[, D] == 0, X]), Y = as.factor(unlist(dat[dat[, D] == 0, S])), seed = seed)
pf1 <- probability_forest(X = as.matrix(dat[dat[, D] == 1, X]), Y = as.factor(unlist(dat[dat[, D] == 1, S])), seed = seed)

pf0_pred <- predict(pf0, newdata = as.matrix(dat[, X]))
pf1_pred <- predict(pf1, newdata = as.matrix(dat[, X]))

one_q0 <- pf0_pred$predictions[, 2]
one_q1 <- pf1_pred$predictions[, 2]

# conditional bounds (not reported)
X_avg <- rep(NA, 26)
X_avg[c(2, 10:13, 23:26)] <- apply(dat[, X[c(2, 10:13, 23:26)]], 2, function(x) mean(x, na.rm = 1))
X_avg[c(1, 3:9, 14:22)] <- 0
Xm_evals <- unique(quantile(dat[, X[27]], seq(0.1, 0.9, 0.1)))
X_moderator <- matrix(rep(X_avg, length(Xm_evals)), length(X_avg), length(Xm_evals))
X_moderator <- rbind(X_moderator, Xm_evals)
X_moderator <- t(X_moderator)

pf0_pred_m <- predict(pf0, newdata = as.matrix(X_moderator))
pf1_pred_m <- predict(pf1, newdata = as.matrix(X_moderator))

one_q0_m <- pf0_pred_m$predictions[, 2]
one_q1_m <- pf1_pred_m$predictions[, 2]
one_q <- one_q0_m / one_q1_m
one_q[one_q > 1] <- 1/one_q[one_q > 1]

result <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S, X = X, 
               W = NULL, Pscore = "Pscore",
               regression.splitting = FALSE, cv_fold = 5,
               aggBounds = 1,
               cBounds = 1, X_moderator = X_moderator, 
               LeeBounds = 1, direction = NULL,
               cond.mono = 0)

tau_l_est_avg <- result[["tau_l_est_avg"]]
tau_u_est_avg <- result[["tau_u_est_avg"]]
se_tau_l_est_avg <- result[["se_tau_l_est_avg"]]
se_tau_u_est_avg <- result[["se_tau_u_est_avg"]]

tau_l_est_lee_mono <- result[["tau_l_est_lee"]]
tau_u_est_lee_mono <- result[["tau_u_est_lee"]]
se_tau_l_est_lee_mono <- result[["se_tau_l_est_lee"]]
se_tau_u_est_lee_mono <- result[["se_tau_u_est_lee"]]

tau_l_m_avg <- result[["tau_l_m_avg"]]
tau_u_m_avg <- result[["tau_u_m_avg"]]
se_tau_l_m_avg <- result[["se_tau_l_m_avg"]]
se_tau_u_m_avg <- result[["se_tau_u_m_avg"]]

result.lee <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S, X = X, 
                   W = NULL, Pscore = "Ps",
                   regression.splitting = FALSE, cv_fold = 5,
                   aggBounds = 0,
                   cBounds = 0, X_moderator = X_moderator, 
                   LeeBounds = 1, direction = NULL,
                   cond.mono = 1)

tau_l_est_lee_condmono <- result.lee[["tau_l_est_lee"]]
tau_u_est_lee_condmono <- result.lee[["tau_u_est_lee"]]
se_tau_l_est_lee_condmono <- result.lee[["se_tau_l_est_lee"]]
se_tau_u_est_lee_condmono <- result.lee[["se_tau_u_est_lee"]]

cf_m <- causal_forest(X = as.matrix(dat[dat[, S] == 1, X]),
                      Y = dat[dat[, S] == 1, Y],
                      W = dat[dat[, S] == 1, D], seed = seed)
pred_m <- predict(cf_m, newdata = as.matrix(X_moderator), estimate.variance = 1)
tau_m_avg <- pred_m[, 1]

pdata_m <- data.frame("X" = X_moderator[, 27], "tau_u_m_avg" = tau_u_m_avg,
                      "tau_l_m_avg" = tau_l_m_avg, "se_tau_u_m_avg" = se_tau_u_m_avg,
                      "se_tau_l_m_avg" = se_tau_l_m_avg)
pdata_m <- pdata_m[order(pdata_m$X), ]
save(pdata_m, file = "application/data/JobCorps_results_m_90.RData")

load("application/data/JobCorps_results_m_90.RData")

pdf("application/graphs/JobCorps_pre_90.pdf", width = 9, height = 9)
par(mar = c(5, 5, 2, 2))
plot(tau_u_m_avg ~ X, pdata_m, type = "p", ylab = "Conditional bounds", 
     xlab = "Pre-treatment wage", col = "red", ylim = c(-.1, .15), pch = 20,
     cex = 2, cex.lab=2, cex.axis=2)
lines(tau_u_m_avg - 1.96*se_tau_u_m_avg ~ X, pdata_m, col = "red", lty = 3, lwd = 3)
lines(tau_u_m_avg + 1.96*se_tau_u_m_avg ~ X, pdata_m, col = "red", lty = 3, lwd = 3)
points(tau_l_m_avg ~ X, pdata_m, col = "blue", lty = 2, pch = 20, cex = 2)
lines(tau_l_m_avg - 1.96*se_tau_l_m_avg ~ X, pdata_m, col = "blue", lty = 3, lwd = 3)
lines(tau_l_m_avg + 1.96*se_tau_l_m_avg ~ X, pdata_m, col = "blue", lty = 3, lwd = 3)
abline(h = 0, lty = 2, col = "gray", lwd = 6)
legend("topleft", pch = c(20, 20, NA, NA), lty = c(NA, NA, 3, 3), lwd = c(NA, NA, 3, 3),
       col = c("red", "blue", "red", "blue"),
       legend = c("Upper bounds", "Lower bounds", "95% CI of upper bounds", "95% CI of lower bounds"), 
       cex = 2, bty = "n")
dev.off()

reg_formula <- as.formula(paste0(Y, "~", D, "+", paste0(X, collapse = "+")))
tau_reg <- lm(reg_formula, dat[dat[, S] == 1, ])
tau_ate <- coef(tau_reg)[2]
se_tau_ate <- sqrt(vcovHC(tau_reg, type = "HC1")[2, 2])
sqrt(vcovCL(tau_reg, cluster = dat$allsides_link)[2, 2])


all_ests <- c(tau_l_est_avg, tau_u_est_avg, 
              tau_l_est_lee_condmono, tau_u_est_lee_condmono,
              tau_l_est_lee_mono, tau_u_est_lee_mono,
              tau_ate)
all_CI95_upper <- c(tau_l_est_avg + 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg + 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono + 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono + 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono + 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono + 1.96*se_tau_u_est_lee_mono,
                    tau_ate + 1.96*se_tau_ate)
all_CI95_lower <- c(tau_l_est_avg - 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg - 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono - 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono - 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono - 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono - 1.96*se_tau_u_est_lee_mono,
                    tau_ate - 1.96*se_tau_ate)


method <- c("CTB", "CTB", "TB (cond. mono.)", "TB (cond. mono.)", "TB (mono.)", "TB (mono.)", "OLS")
method <- factor(method, levels = c("OLS", "TB (mono.)", "TB (cond. mono.)", "CTB"))

type <- c("Lower", "Upper", "Lower", "Upper", "Lower", "Upper", "ATE")

p.data <- data.frame("Estimates" = round(all_ests, 3), "CI95_u" = round(all_CI95_upper, 3), 
                     "CI95_l" = round(all_CI95_lower, 3), "Method" = method, "Estimand" = type)
save(p.data, file = "application/data/JobCorps_results_agg_90.RData")

load("application/data/JobCorps_results_agg_90.RData")

pd <- position_dodge(width=0.2)

coefs_p <- ggplot(p.data, aes(Method, Estimates, color=Estimand)) +
  geom_point(aes(shape=Method), size=4, position=pd, show.legend = FALSE) +
  scale_color_manual(name="Estimand",values=c("black", "blue", "red")) +
  scale_shape_manual(name="Method",values=c(19, 17, 17, 15)) +
  theme_bw() +
  geom_errorbar(aes(Method, ymin=CI95_l, ymax=CI95_u), width=0.1, position=pd) + 
  geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('') + ylab('') +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), axis.text.x = element_text(angle = 0), legend.position="top") + 
  ylab("Impact on wage in week 90") +
  coord_flip()
ggsave(paste0("application/graphs/compareBounds_JobCorps_90.pdf"), coefs_p, device = "pdf", height = 9, width = 9) 


# week 208
leedata=data.frame(treat=JobCorps_data_baseline$TREATMNT.y,
                   selection=JobCorps_data_employment$week_208,
                   outcome=JobCorps_data_wages$week_208,
                   JobCorps_data_baseline[, 3:29])
leedata$logwage <- log(leedata$WKEARNR+1)

Y <- "outcome"
D <- "treat"
S <- "selection"
X <- names(leedata)[c(4:29, 31)]
Pscore <- "Pscore"

dat <- leedata[complete.cases(leedata[, D]), ]
dat <- dat[dat[, Y] != Inf, ]
for (i in 1:length(X)){
  Xi <- X[i]
  dat[, Xi] <- as.numeric(dat[, Xi])
}

Ntr <- sum(dat[, D] == 1)
N <- nrow(dat)

dat$Pscore <- Ntr/N

pf0 <- probability_forest(X = as.matrix(dat[dat[, D] == 0, X]), Y = as.factor(unlist(dat[dat[, D] == 0, S])), seed = seed)
pf1 <- probability_forest(X = as.matrix(dat[dat[, D] == 1, X]), Y = as.factor(unlist(dat[dat[, D] == 1, S])), seed = seed)

pf0_pred <- predict(pf0, newdata = as.matrix(dat[, X]))
pf1_pred <- predict(pf1, newdata = as.matrix(dat[, X]))

one_q0 <- pf0_pred$predictions[, 2]
one_q1 <- pf1_pred$predictions[, 2]

# conditional bounds (not reported)
X_avg <- rep(NA, 26)
X_avg[c(2, 10:13, 23:26)] <- apply(dat[, X[c(2, 10:13, 23:26)]], 2, function(x) mean(x, na.rm = 1))
X_avg[c(1, 3:9, 14:22)] <- 0
Xm_evals <- unique(quantile(dat[, X[27]], seq(0.1, 0.9, 0.1)))
X_moderator <- matrix(rep(X_avg, length(Xm_evals)), length(X_avg), length(Xm_evals))
X_moderator <- rbind(X_moderator, Xm_evals)
X_moderator <- t(X_moderator)

pf0_pred_m <- predict(pf0, newdata = as.matrix(X_moderator))
pf1_pred_m <- predict(pf1, newdata = as.matrix(X_moderator))

one_q0_m <- pf0_pred_m$predictions[, 2]
one_q1_m <- pf1_pred_m$predictions[, 2]
one_q <- one_q0_m / one_q1_m
one_q[one_q > 1] <- 1/one_q[one_q > 1]

result <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S, X = X, 
               W = NULL, Pscore = "Pscore",
               regression.splitting = FALSE, cv_fold = 5,
               aggBounds = 1,
               cBounds = 1, X_moderator = X_moderator, 
               LeeBounds = 1, direction = NULL,
               cond.mono = 0)

tau_l_est_avg <- result[["tau_l_est_avg"]]
tau_u_est_avg <- result[["tau_u_est_avg"]]
se_tau_l_est_avg <- result[["se_tau_l_est_avg"]]
se_tau_u_est_avg <- result[["se_tau_u_est_avg"]]

tau_l_est_lee_mono <- result[["tau_l_est_lee"]]
tau_u_est_lee_mono <- result[["tau_u_est_lee"]]
se_tau_l_est_lee_mono <- result[["se_tau_l_est_lee"]]
se_tau_u_est_lee_mono <- result[["se_tau_u_est_lee"]]

tau_l_m_avg <- result[["tau_l_m_avg"]]
tau_u_m_avg <- result[["tau_u_m_avg"]]
se_tau_l_m_avg <- result[["se_tau_l_m_avg"]]
se_tau_u_m_avg <- result[["se_tau_u_m_avg"]]

result.lee <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S, X = X, 
                   W = NULL, Pscore = "Ps",
                   regression.splitting = FALSE, cv_fold = 5,
                   aggBounds = 0,
                   cBounds = 0, X_moderator = X_moderator, 
                   LeeBounds = 1, direction = NULL,
                   cond.mono = 1)

tau_l_est_lee_condmono <- result.lee[["tau_l_est_lee"]]
tau_u_est_lee_condmono <- result.lee[["tau_u_est_lee"]]
se_tau_l_est_lee_condmono <- result.lee[["se_tau_l_est_lee"]]
se_tau_u_est_lee_condmono <- result.lee[["se_tau_u_est_lee"]]

cf_m <- causal_forest(X = as.matrix(dat[dat[, S] == 1, X]),
                      Y = dat[dat[, S] == 1, Y],
                      W = dat[dat[, S] == 1, D], seed = seed)
pred_m <- predict(cf_m, newdata = as.matrix(X_moderator), estimate.variance = 1)
tau_m_avg <- pred_m[, 1]

pdata_m <- data.frame("X" = X_moderator[, 27], "tau_u_m_avg" = tau_u_m_avg,
                      "tau_l_m_avg" = tau_l_m_avg, "se_tau_u_m_avg" = se_tau_u_m_avg,
                      "se_tau_l_m_avg" = se_tau_l_m_avg)
pdata_m <- pdata_m[order(pdata_m$X), ]
save(pdata_m, file = "application/data/JobCorps_results_m_208.RData")

load("application/data/JobCorps_results_m_208.RData")

pdf("application/graphs/JobCorps_pre_208.pdf", width = 9, height = 9)
par(mar = c(5, 5, 2, 2))
plot(tau_u_m_avg ~ X, pdata_m, type = "p", ylab = "Conditional bounds", 
     xlab = "Pre-treatment wage", col = "red", ylim = c(-.1, .15), pch = 20,
     cex = 2, cex.lab=2, cex.axis=2)
lines(tau_u_m_avg - 1.96*se_tau_u_m_avg ~ X, pdata_m, col = "red", lty = 3, lwd = 3)
lines(tau_u_m_avg + 1.96*se_tau_u_m_avg ~ X, pdata_m, col = "red", lty = 3, lwd = 3)
points(tau_l_m_avg ~ X, pdata_m, col = "blue", lty = 2, pch = 20, cex = 2)
lines(tau_l_m_avg - 1.96*se_tau_l_m_avg ~ X, pdata_m, col = "blue", lty = 3, lwd = 3)
lines(tau_l_m_avg + 1.96*se_tau_l_m_avg ~ X, pdata_m, col = "blue", lty = 3, lwd = 3)
abline(h = 0, lty = 2, col = "gray", lwd = 6)
legend("topleft", pch = c(20, 20, NA, NA), lty = c(NA, NA, 3, 3), lwd = c(NA, NA, 3, 3),
       col = c("red", "blue", "red", "blue"),
       legend = c("Upper bounds", "Lower bounds", "95% CI of upper bounds", "95% CI of lower bounds"), 
       cex = 2, bty = "n")
dev.off()

reg_formula <- as.formula(paste0(Y, "~", D, "+", paste0(X, collapse = "+")))
tau_reg <- lm(reg_formula, dat[dat[, S] == 1, ])
tau_ate <- coef(tau_reg)[2]
se_tau_ate <- sqrt(vcovHC(tau_reg, type = "HC1")[2, 2])
sqrt(vcovCL(tau_reg, cluster = dat$allsides_link)[2, 2])


all_ests <- c(tau_l_est_avg, tau_u_est_avg, 
              tau_l_est_lee_condmono, tau_u_est_lee_condmono,
              tau_l_est_lee_mono, tau_u_est_lee_mono,
              tau_ate)
all_CI95_upper <- c(tau_l_est_avg + 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg + 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono + 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono + 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono + 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono + 1.96*se_tau_u_est_lee_mono,
                    tau_ate + 1.96*se_tau_ate)
all_CI95_lower <- c(tau_l_est_avg - 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg - 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono - 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono - 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono - 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono - 1.96*se_tau_u_est_lee_mono,
                    tau_ate - 1.96*se_tau_ate)


method <- c("CTB", "CTB", "TB (cond. mono.)", "TB (cond. mono.)", "TB (mono.)", "TB (mono.)", "OLS")
method <- factor(method, levels = c("OLS", "TB (mono.)", "TB (cond. mono.)", "CTB"))

type <- c("Lower", "Upper", "Lower", "Upper", "Lower", "Upper", "ATE")

p.data <- data.frame("Estimates" = round(all_ests, 3), "CI95_u" = round(all_CI95_upper, 3), 
                     "CI95_l" = round(all_CI95_lower, 3), "Method" = method, "Estimand" = type)
save(p.data, file = "application/data/JobCorps_results_agg_208.RData")

load("application/data/JobCorps_results_agg_208.RData")

pd <- position_dodge(width=0.2)

coefs_p <- ggplot(p.data, aes(Method, Estimates, color=Estimand)) +
  geom_point(aes(shape=Method), size=4, position=pd, show.legend = FALSE) +
  scale_color_manual(name="Estimand",values=c("black", "blue", "red")) +
  scale_shape_manual(name="Method",values=c(19, 17, 17, 15)) +
  theme_bw() +
  geom_errorbar(aes(Method, ymin=CI95_l, ymax=CI95_u), width=0.1, position=pd) + 
  geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('') + ylab('') +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), axis.text.x = element_text(angle = 0), legend.position="top") + 
  ylab("Impact on wage in week 90") +
  coord_flip()
ggsave(paste0("application/graphs/compareBounds_JobCorps_208.pdf"), coefs_p, device = "pdf", height = 9, width = 9) 

