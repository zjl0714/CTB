######################################################
##### Replication of Blattman and Annan (2010) #######
######################################################

rm(list=ls())
library(dplyr)
library(grf)
library(hdm)
library(ggplot2)
library(sandwich)
library(foreign)

setwd("/Users/yewang/Documents/GitHub/trimming/replication/")
source("function/CTTB.R")

seed <- 42

dat <- read.dta("application/data/SWAY_I_June2013.dta")

AC <- names(dat)[c(280:286, 288:303)] # age dummies and race dummies
ACG <- c(AC, "G1_18_21", "G1_22_25", "G2_14_17", "G2_22_25", "G3_14_17", "G3_18_21", "G4_14_17", "G4_18_21", "G4_22_25", 
         "G5_14_17", "G5_18_21", "G5_22_25", "G6_14_17", "G6_18_21", "G6_22_25", "G7_14_17", "G7_18_21", "G7_22_25", "G8_14_17", "G8_18_21", "G8_22_25")
X <- c("hh_fthr_frm", "hh_size96", "hh_land", "hh_cattle", "hh_stock", "hh_plow", "age", ACG)

# education and sampling attrition
Y <- "educ"
S1 <- "sample"
S2 <- "found"
D <- "abd"
W <- "i_weight"
W1 <- "w_abs2"
W2 <- "w_sel_abs"

dat <- dat[complete.cases(dat[, c(D, S1, X, W)]), ]
dat <- dat[, c(Y, D, S1, S2, X, W, W1, W2)]

# estimate q0 and q1
pf0 <- probability_forest(X = as.matrix(dat[dat[, D] == 0, X]), Y = as.factor(unlist(dat[dat[, D] == 0, S1])), seed = seed)
pf1 <- probability_forest(X = as.matrix(dat[dat[, D] == 1, X]), Y = as.factor(unlist(dat[dat[, D] == 1, S1])), seed = seed)

pf0_pred <- predict(pf0, newdata = as.matrix(dat[, X]))
pf1_pred <- predict(pf1, newdata = as.matrix(dat[, X]))

one_q0 <- pf0_pred$predictions[, 2]
one_q1 <- pf1_pred$predictions[, 2]

pdf("application/graphs/cBounds_BA_2010_qs1.pdf", width = 9, height = 9)
par(mar = c(5, 5, 2, 2))
plot(density(one_q0), ylab = "Density", main = expression("S"[1]),
     xlab = "Response rate", col = "blue", ylim = c(0, 8), xlim = c(0, 1),
     cex = 2, cex.lab=2, cex.axis=2, cex.main=2)
lines(density(one_q1), col = "red", lty = 1)
legend("topleft", lty = c(1, 1), 
       col = c("red", "blue", "black"),
       legend = c("Response rate (D = 1)", "Response rate (D = 0)"), 
       cex = 2, bty = "n")
dev.off()

pdf("application/graphs/cBounds_BA_2010_qratio1.pdf", width = 9, height = 9)
par(mar = c(5, 5, 2, 2))
plot(density(one_q0 / one_q1), ylab = "Density", main = expression("S"[1]),
     xlab = "Relative share of always-responders", col = "black", ylim = c(0, 4), xlim = c(0, 2),
     cex = 2, cex.lab=2, cex.axis=2, cex.main=2)
abline(v = 1, lty = 3, col = "gray", lwd = 6)
dev.off()

pf0 <- probability_forest(X = as.matrix(dat[dat[, D] == 0, X]), Y = as.factor(unlist(dat[dat[, D] == 0, S2])), seed = seed)
pf1 <- probability_forest(X = as.matrix(dat[dat[, D] == 1, X]), Y = as.factor(unlist(dat[dat[, D] == 1, S2])), seed = seed)

pf0_pred <- predict(pf0, newdata = as.matrix(dat[, X]))
pf1_pred <- predict(pf1, newdata = as.matrix(dat[, X]))

one_q0 <- pf0_pred$predictions[, 2]
one_q1 <- pf1_pred$predictions[, 2]

pdf("application/graphs/cBounds_BA_2010_qs2.pdf", width = 9, height = 9)
par(mar = c(5, 5, 2, 2))
plot(density(one_q0), ylab = "Density", main = expression("S"[2]),
     xlab = "Response rate", col = "blue", ylim = c(0, 8), xlim = c(0, 1),
     cex = 2, cex.lab=2, cex.axis=2, cex.main=2)
lines(density(one_q1), col = "red", lty = 1)
legend("topleft", lty = c(1, 1), 
       col = c("red", "blue", "black"),
       legend = c("Response rate (D = 1)", "Response rate (D = 0)"), 
       cex = 2, bty = "n")
dev.off()

pdf("application/graphs/cBounds_BA_2010_qratio2.pdf", width = 9, height = 9)
par(mar = c(5, 5, 2, 2))
plot(density(one_q0 / one_q1), ylab = "Density", main = expression("S"[2]),
     xlab = "Relative share of always-responders", col = "black", ylim = c(0, 4), xlim = c(0, 2),
     cex = 2, cex.lab=2, cex.axis=2, cex.main=2)
abline(v = 1, lty = 3, col = "gray", lwd = 6)
dev.off()

result <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S1, X = X, W = W, Pscore = NULL,
               regression.splitting = FALSE, cv_fold = 10,
               aggBounds = 1, trim_l = 0.05, trim_u = 0.95,
               cBounds = 0, X_moderator = NULL, 
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

result <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S1, X = X, W = W, Pscore = NULL,
               regression.splitting = FALSE, cv_fold = 10,
               aggBounds = 0,
               cBounds = 0, X_moderator = NULL, 
               LeeBounds = 1, direction = NULL,
               cond.mono = 1)

tau_l_est_lee_condmono <- result[["tau_l_est_lee"]]
tau_u_est_lee_condmono <- result[["tau_u_est_lee"]]
se_tau_l_est_lee_condmono <- result[["se_tau_l_est_lee"]]
se_tau_u_est_lee_condmono <- result[["se_tau_u_est_lee"]]

# estimate the ate using observations that select into the sample
## regression
ate_reg <- lm(as.formula(paste0(Y, "~", D, "+", paste0(X, collapse = "+"))), data = dat, weights = w_sel_abs)
tau_ate <- coef(ate_reg)[2]
se_tau_ate <- summary(ate_reg)$coefficients[2, 2]

## causal forest with ipw
ate_cf <- causal_forest(X = as.matrix(dat[dat[, S2] == 1, X]), 
                        Y = c(dat[dat[, S2] == 1, Y]), 
                        W = c(dat[dat[, S2] == 1, D]), 
                        sample.weights = c(dat[dat[, S2] == 1, W1]), seed = seed)
ate_est_cf <- average_treatment_effect(ate_cf)
tau_ate_cf <- ate_est_cf[1]
se_tau_ate_cf <- ate_est_cf[2]

## causal forest without ipw
ate_cf_naive <- causal_forest(X = as.matrix(dat[dat[, S2] == 1, X]), 
                              Y = c(dat[dat[, S2] == 1, Y]), 
                              W = c(dat[dat[, S2] == 1, D]))
ate_est_cf_naive <- average_treatment_effect(ate_cf_naive)
tau_ate_cf_naive <- ate_est_cf_naive[1]
se_tau_ate_cf_naive <- ate_est_cf_naive[2]
  
all_ests <- c(tau_l_est_avg, tau_u_est_avg, 
              tau_l_est_lee_condmono, tau_u_est_lee_condmono,
              tau_l_est_lee_mono, tau_u_est_lee_mono,
              tau_ate_cf, tau_ate_cf_naive)
all_CI95_upper <- c(tau_l_est_avg + 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg + 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono + 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono + 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono + 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono + 1.96*se_tau_u_est_lee_mono,
                    tau_ate_cf + 1.96*se_tau_ate_cf,
                    tau_ate_cf_naive + 1.96*se_tau_ate_cf_naive)
all_CI95_lower <- c(tau_l_est_avg - 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg - 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono - 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono - 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono - 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono - 1.96*se_tau_u_est_lee_mono,
                    tau_ate_cf - 1.96*se_tau_ate_cf,
                    tau_ate_cf_naive - 1.96*se_tau_ate_cf_naive)


method <- c("CTB", "CTB", "TB (cond. mono.)", "TB (cond. mono.)", "TB (mono.)", "TB (mono.)", "CF", "CF (no correction)")
method <- factor(method, levels = c("CF (no correction)", "CF", "TB (mono.)", "TB (cond. mono.)", "CTB"))

type <- c("Lower", "Upper", "Lower", "Upper", "Lower", "Upper", "ATE", "ATE")

p.data <- data.frame("Estimates" = round(all_ests, 3), "CI95_u" = round(all_CI95_upper, 3), 
                     "CI95_l" = round(all_CI95_lower, 3), "Method" = method, "Estimand" = type)
save(p.data, file = "application/data/Blattman_Annan_2010_results_educ.RData")

load("application/data/Blattman_Annan_2010_results_educ.RData")

pd <- position_dodge(width=0.2)

coefs_p <- ggplot(p.data, aes(Method, Estimates, color=Estimand)) +
  geom_point(aes(shape=Method), size=4, position=pd, show.legend = FALSE) +
  scale_color_manual(name="Estimand",values=c("black", "blue", "red")) +
  scale_shape_manual(name="Method",values=c(19, 19, 17, 17, 15)) +
  theme_bw() +
  geom_errorbar(aes(Method, ymin=CI95_l, ymax=CI95_u), width=0.1, position=pd) + 
  geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('') + ylab('') +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), axis.text.x = element_text(angle = 0), legend.position = "top") + 
  ylab("Impact on years of education") +
  coord_flip()
ggsave(paste0("application/graphs/compareBounds_BA_2010_educ.pdf"), coefs_p, device = "pdf", height = 9, width = 9) 



# distress
dat <- read.dta("application/data/SWAY_I_June2013.dta")

AC <- names(dat)[c(280:286, 288:303)] # age dummies and race dummies
ACG <- c(AC, "G1_18_21", "G1_22_25", "G2_14_17", "G2_22_25", "G3_14_17", "G3_18_21", "G4_14_17", "G4_18_21", "G4_22_25", 
         "G5_14_17", "G5_18_21", "G5_22_25", "G6_14_17", "G6_18_21", "G6_22_25", "G7_14_17", "G7_18_21", "G7_22_25", "G8_14_17", "G8_18_21", "G8_22_25")
X <- c("hh_fthr_frm", "hh_size96", "hh_land", "hh_cattle", "hh_stock", "hh_plow", "age", ACG)

Y <- "distress"
S <- "found"
D <- "abd"
W <- "i_weight"
W1 <- "w_abs2"
W2 <- "w_sel_abs"

dat <- dat[complete.cases(dat[, c(D, S, X, W)]), ]
dat <- dat[, c(Y, D, S, X, W, W1, W2)]

result <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S, X = X, W = W, Pscore = NULL,
               regression.splitting = FALSE, cv_fold = 10,
               aggBounds = 1, trim_l = 0.05, trim_u = 0.95,
               cBounds = 0, X_moderator = NULL, 
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

result <- CTTB(data = dat, Y = Y, D = D, S = S, X = X, W = W, Pscore = NULL,
               regression.splitting = FALSE, cv_fold = 10,
               aggBounds = 0,
               cBounds = 0, X_moderator = NULL, 
               LeeBounds = 1, direction = NULL,
               cond.mono = 1)

tau_l_est_lee_condmono <- result[["tau_l_est_lee"]]
tau_u_est_lee_condmono <- result[["tau_u_est_lee"]]
se_tau_l_est_lee_condmono <- result[["se_tau_l_est_lee"]]
se_tau_u_est_lee_condmono <- result[["se_tau_u_est_lee"]]

# estimate the ate using observations that select into the sample
## regression
ate_reg <- lm(as.formula(paste0(Y, "~", D, "+", paste0(X, collapse = "+"))), data = dat, weights = i_weight)
tau_ate <- coef(ate_reg)[2]
se_tau_ate <- summary(ate_reg)$coefficients[2, 2]

## causal forest with ipw
ate_cf <- causal_forest(X = as.matrix(dat[dat[, S] == 1, X]), 
                        Y = c(dat[dat[, S] == 1, Y]), 
                        W = c(dat[dat[, S] == 1, D]), 
                        sample.weights = c(dat[dat[, S] == 1, W1]),
                        seed = seed)
ate_est_cf <- average_treatment_effect(ate_cf)
tau_ate_cf <- ate_est_cf[1]
se_tau_ate_cf <- ate_est_cf[2]

## causal forest without ipw
ate_cf_naive <- causal_forest(X = as.matrix(dat[dat[, S] == 1, X]), 
                              Y = c(dat[dat[, S] == 1, Y]), 
                              W = c(dat[dat[, S] == 1, D]),
                              seed = seed)
ate_est_cf_naive <- average_treatment_effect(ate_cf_naive)
tau_ate_cf_naive <- ate_est_cf_naive[1]
se_tau_ate_cf_naive <- ate_est_cf_naive[2]

all_ests <- c(tau_l_est_avg, tau_u_est_avg, 
              tau_l_est_lee_condmono, tau_u_est_lee_condmono,
              tau_l_est_lee_mono, tau_u_est_lee_mono,
              tau_ate_cf, tau_ate_cf_naive)
all_CI95_upper <- c(tau_l_est_avg + 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg + 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono + 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono + 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono + 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono + 1.96*se_tau_u_est_lee_mono,
                    tau_ate_cf + 1.96*se_tau_ate_cf,
                    tau_ate_cf_naive + 1.96*se_tau_ate_cf_naive)
all_CI95_lower <- c(tau_l_est_avg - 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg - 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono - 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono - 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono - 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono - 1.96*se_tau_u_est_lee_mono,
                    tau_ate_cf - 1.96*se_tau_ate_cf,
                    tau_ate_cf_naive - 1.96*se_tau_ate_cf_naive)


method <- c("CTB", "CTB", "TB (cond. mono.)", "TB (cond. mono.)", "TB (mono.)", "TB (mono.)", "CF", "CF (no correction)")
method <- factor(method, levels = c("CF (no correction)", "CF", "TB (mono.)", "TB (cond. mono.)", "CTB"))

type <- c("Lower", "Upper", "Lower", "Upper", "Lower", "Upper", "ATE", "ATE")

p.data <- data.frame("Estimates" = round(all_ests, 3), "CI95_u" = round(all_CI95_upper, 3), 
                     "CI95_l" = round(all_CI95_lower, 3), "Method" = method, "Estimand" = type)
save(p.data, file = "application/data/Blattman_Annan_2010_results_distress.RData")

load("application/data/Blattman_Annan_2010_results_distress.RData")

pd <- position_dodge(width=0.2)

coefs_p <- ggplot(p.data, aes(Method, Estimates, color=Estimand)) +
  geom_point(aes(shape=Method), size=4, position=pd, show.legend = FALSE) +
  scale_color_manual(name="Estimand",values=c("black", "blue", "red")) +
  scale_shape_manual(name="Method",values=c(19, 19, 17, 17, 15)) +
  theme_bw() +
  geom_errorbar(aes(Method, ymin=CI95_l, ymax=CI95_u), width=0.1, position=pd) + 
  geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('') + ylab('') +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), axis.text.x = element_text(angle = 0), legend.position = "top") + 
  ylab("Impact on distress") +
  coord_flip()
ggsave(paste0("application/graphs/compareBounds_BA_2010_distress.pdf"), coefs_p, device = "pdf", height = 9, width = 9) 

# wage
dat <- read.dta("application/data/SWAY_I_June2013.dta")

AC <- names(dat)[c(280:286, 288:303)] # age dummies and race dummies
ACG <- c(AC, "G1_18_21", "G1_22_25", "G2_14_17", "G2_22_25", "G3_14_17", "G3_18_21", "G4_14_17", "G4_18_21", "G4_22_25", 
         "G5_14_17", "G5_18_21", "G5_22_25", "G6_14_17", "G6_18_21", "G6_22_25", "G7_14_17", "G7_18_21", "G7_22_25", "G8_14_17", "G8_18_21", "G8_22_25")
X <- c("hh_fthr_frm", "hh_size96", "hh_land", "hh_cattle", "hh_stock", "hh_plow", "age", ACG)

dat$has_wage <- !is.na(dat$wage_mo)

Y <- "wage_mo"
S <- "has_wage"
D <- "abd"
W <- "i_weight"
W1 <- "w_abs2"
W2 <- "w_sel_abs"

dat <- dat[complete.cases(dat[, c(D, S, X, W)]), ]
dat <- dat[, c(Y, D, S, X, W, W1, W2)]


result <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S, X = X, W = W, 
               Pscore = NULL, regression.splitting = FALSE, cv_fold = 10, 
               trim_l = 0.05, trim_u = 0.95, 
               aggBounds = 1,
               cBounds = 0, X_moderator = NULL, 
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

result <- CTTB(data = dat, Y = Y, D = D, S = S, X = X, W = W, Pscore = NULL,
               regression.splitting = FALSE, cv_fold = 10,
               aggBounds = 0,
               cBounds = 0, X_moderator = NULL, 
               LeeBounds = 1, direction = NULL,
               cond.mono = 1)

tau_l_est_lee_condmono <- result[["tau_l_est_lee"]]
tau_u_est_lee_condmono <- result[["tau_u_est_lee"]]
se_tau_l_est_lee_condmono <- result[["se_tau_l_est_lee"]]
se_tau_u_est_lee_condmono <- result[["se_tau_u_est_lee"]]

ate_reg <- lm(as.formula(paste0(Y, "~", D, "+", paste0(X, collapse = "+"))), data = dat, weights = i_weight)
tau_ate <- coef(ate_reg)[2]
se_tau_ate <- summary(ate_reg)$coefficients[2, 2]

ate_cf <- causal_forest(X = as.matrix(dat[dat[, S] == 1, X]), 
                        Y = c(dat[dat[, S] == 1, Y]), 
                        W = c(dat[dat[, S] == 1, D]), 
                        sample.weights = c(dat[dat[, S] == 1, W1]),
                        seed = seed)
ate_est_cf <- average_treatment_effect(ate_cf)
tau_ate_cf <- ate_est_cf[1]
se_tau_ate_cf <- ate_est_cf[2]

ate_cf_naive <- causal_forest(X = as.matrix(dat[dat[, S] == 1, X]), 
                              Y = c(dat[dat[, S] == 1, Y]), 
                              W = c(dat[dat[, S] == 1, D]))
ate_est_cf_naive <- average_treatment_effect(ate_cf_naive)
tau_ate_cf_naive <- ate_est_cf_naive[1]
se_tau_ate_cf_naive <- ate_est_cf_naive[2]

all_ests <- c(tau_l_est_avg, tau_u_est_avg, 
              tau_l_est_lee_condmono, tau_u_est_lee_condmono,
              tau_l_est_lee_mono, tau_u_est_lee_mono,
              tau_ate_cf, tau_ate_cf_naive)
all_CI95_upper <- c(tau_l_est_avg + 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg + 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono + 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono + 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono + 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono + 1.96*se_tau_u_est_lee_mono,
                    tau_ate_cf + 1.96*se_tau_ate_cf,
                    tau_ate_cf_naive + 1.96*se_tau_ate_cf_naive)
all_CI95_lower <- c(tau_l_est_avg - 1.96*se_tau_l_est_avg, 
                    tau_u_est_avg - 1.96*se_tau_u_est_avg, 
                    tau_l_est_lee_condmono - 1.96*se_tau_l_est_lee_condmono, 
                    tau_u_est_lee_condmono - 1.96*se_tau_u_est_lee_condmono,
                    tau_l_est_lee_mono - 1.96*se_tau_l_est_lee_mono, 
                    tau_u_est_lee_mono - 1.96*se_tau_u_est_lee_mono,
                    tau_ate_cf - 1.96*se_tau_ate_cf,
                    tau_ate_cf_naive - 1.96*se_tau_ate_cf_naive)


method <- c("CTB", "CTB", "TB (cond. mono.)", "TB (cond. mono.)", "TB (mono.)", "TB (mono.)", "CF", "CF (no correction)")
method <- factor(method, levels = c("CF (no correction)", "CF", "TB (mono.)", "TB (cond. mono.)", "CTB"))

type <- c("Lower", "Upper", "Lower", "Upper", "Lower", "Upper", "ATE", "ATE")

p.data <- data.frame("Estimates" = round(all_ests, 3), "CI95_u" = round(all_CI95_upper, 3), 
                     "CI95_l" = round(all_CI95_lower, 3), "Method" = method, "Estimand" = type)
save(p.data, file = "application/data/Blattman_Annan_2010_results_wage.RData")

load("application/data/Blattman_Annan_2010_results_wage.RData")

pd <- position_dodge(width=0.2)

coefs_p <- ggplot(p.data, aes(Method, Estimates, color=Estimand)) +
  geom_point(aes(shape=Method), size=4, position=pd, show.legend = FALSE) +
  scale_color_manual(name="Estimand",values=c("black", "blue", "red")) +
  scale_shape_manual(name="Method",values=c(19, 19, 17, 17, 15)) +
  theme_bw() +
  geom_errorbar(aes(Method, ymin=CI95_l, ymax=CI95_u), width=0.1, position=pd) + 
  geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('') + ylab('') +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), axis.text.x = element_text(angle = 0), legend.position = "top") + 
  ylab("Impact on wage") +
  coord_flip()
ggsave(paste0("application/graphs/compareBounds_BA_2010_wage.pdf"), coefs_p, device = "pdf", height = 9, width = 9) 

