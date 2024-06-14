######################################################
##### Replication of Kalla and Broockman (2022) ######
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

dat <- get(load("application/data/for_analysis.RData"))

t0.covariate <- c("t0_factor_immpolicy", "t0_factor_immprej",
                  "t0_factor_asyl_policy", "t0_factor_asyl_prej",
                  "t0_factor_lgbtpolicy", "t0_factor_lgbtprej",
                  "t0_therm_trans", "t0_therm_gay", "t0_therm_immigrant", "t0_therm_asylum",
                  "baseline_factor_general_politics",
                  "t0_economy_today", "t0_economy_yearfromnow", "t0_economy_personal",
                  "t0_knowledge_imm_taxes", "t0_college_educ", "t0_ideology", "t0_pid7",
                  "t0_trump_2020", "t0_congress_2020",
                  "t0_employ_fulltime", "t0_employ_parttime", "t0_live_city",
                  "vf_age", "vf_general18", "vf_general16", "vf_general14", "vf_general12",
                  "vf_dem", "vf_rep", "vf_female", "vf_white", "vf_afam", "vf_hispanic",
                  "state_CA", "state_CO", "state_MI", "state_NC", "state_TN")

X <- t0.covariate
Y <- "t1_recall_imm"
S <- "t1_respondent"
Ds <- c("treat_eddie", "treat_maricruz")
D <- Ds[2]

dat <- dat[dat[, Ds[1]] == 0, ]

dat$Ps <- mean(dat[, D])
Ps <- "Ps"

dat <- dat[complete.cases(dat[, c(D, S, X, Ps)]), ]


# CTB
result <- CTTB(data = dat, seed = seed, Y = Y, D = D, S = S, X = X, W = NULL, 
               Pscore = "Ps",
               regression.splitting = FALSE, cv_fold = 5,
               aggBounds = 1,
               cBounds = 0, X_moderator = NULL,
               LeeBounds = 0, direction = NULL,
               cond.mono = 0)
tau_l_est_avg <- result[["tau_l_est_avg"]]
tau_u_est_avg <- result[["tau_u_est_avg"]]
se_tau_l_est_avg <- result[["se_tau_l_est_avg"]]
se_tau_u_est_avg <- result[["se_tau_u_est_avg"]]


# Lee bounds
result <- CTTB(data = dat, Y = Y, D = D, S = S, X = X, W = NULL, Pscore = Ps,
               regression.splitting = FALSE, cv_fold = 5,
               aggBounds = 0,
               cBounds = 0, X_moderator = NULL,
               LeeBounds = 1, direction = NULL,
               cond.mono = 1)

tau_l_est_lee_condmono <- result[["tau_l_est_lee"]]
tau_u_est_lee_condmono <- result[["tau_u_est_lee"]]
se_tau_l_est_lee_condmono <- result[["se_tau_l_est_lee"]]
se_tau_u_est_lee_condmono <- result[["se_tau_u_est_lee"]]

result <- CTTB(data = dat, Y = Y, D = D, S = S, X = X, W = NULL, Pscore = NULL,
               regression.splitting = FALSE, cv_fold = 5,
               aggBounds = 0,
               cBounds = 0, X_moderator = NULL,
               LeeBounds = 1, direction = NULL,
               cond.mono = 0)

tau_l_est_lee_mono <- result[["tau_l_est_lee"]]
tau_u_est_lee_mono <- result[["tau_u_est_lee"]]
se_tau_l_est_lee_mono <- result[["se_tau_l_est_lee"]]
se_tau_u_est_lee_mono <- result[["se_tau_u_est_lee"]]

ate_cf_naive <- causal_forest(X = as.matrix(dat[dat[, S] == 1, X]), 
                              Y = c(dat[dat[, S] == 1, Y]), 
                              W = c(dat[dat[, S] == 1, D]),
                              seed = seed)
ate_est_cf_naive <- average_treatment_effect(ate_cf_naive)
tau_ate_cf_naive <- ate_est_cf_naive[1]
se_tau_ate_cf_naive <- ate_est_cf_naive[2]

all_ests <- c(tau_l_est_avg, tau_u_est_avg,
              tau_l_est_lee_mono, tau_u_est_lee_mono,
              tau_l_est_lee_condmono, tau_u_est_lee_condmono,
              tau_ate_cf_naive)
all_CI95_upper <- c(tau_l_est_avg + 1.96*se_tau_l_est_avg,
                    tau_u_est_avg + 1.96*se_tau_u_est_avg,
                    tau_l_est_lee_mono + 1.96*se_tau_l_est_lee_mono,
                    tau_u_est_lee_mono + 1.96*se_tau_u_est_lee_mono,
                    tau_l_est_lee_condmono + 1.96*se_tau_l_est_lee_condmono,
                    tau_u_est_lee_condmono + 1.96*se_tau_u_est_lee_condmono,
                    tau_ate_cf_naive + 1.96*se_tau_ate_cf_naive)
all_CI95_lower <- c(tau_l_est_avg - 1.96*se_tau_l_est_avg,
                    tau_u_est_avg - 1.96*se_tau_u_est_avg,
                    tau_l_est_lee_mono - 1.96*se_tau_l_est_lee_mono,
                    tau_u_est_lee_mono - 1.96*se_tau_u_est_lee_mono,
                    tau_l_est_lee_condmono - 1.96*se_tau_l_est_lee_condmono,
                    tau_u_est_lee_condmono - 1.96*se_tau_u_est_lee_condmono,
                    tau_ate_cf_naive - 1.96*se_tau_ate_cf_naive)



# plot
method <- c("CTB", "CTB", "TB (mono.)", "TB (mono.)", "TB (cond. mono.)", "TB (cond. mono.)", "CF (no correction)")
method <- factor(method, levels = c("CF (no correction)", "TB (cond. mono.)", "TB (mono.)", "CTB"))

type <- c("Lower", "Upper", "Lower", "Upper", "Lower", "Upper", "ATE")

p.data <- data.frame("Estimates" = round(all_ests, 3), "CI95_u" = round(all_CI95_upper, 3), 
                     "CI95_l" = round(all_CI95_lower, 3), "Method" = method, "Estimand" = type)
save(p.data, file = "application/data/Kalla_Broockman_2022_results_agg.RData")

load("application/data/Kalla_Broockman_2022_results_agg.RData")

pd <- position_dodge(width=0.2)

coefs_p <- ggplot(p.data, aes(Method, Estimates, color=Estimand)) +
  geom_point(aes(shape=Method), size=4, position=pd, show.legend = FALSE) +
  scale_color_manual(name="Estimand", values=c("black", "blue", "red")) +
  scale_shape_manual(name="Method", values=c(19, 19, 17, 17, 15)) +
  theme_bw() +
  geom_errorbar(aes(Method, ymin=CI95_l, ymax=CI95_u), width=0.1, position=pd) + 
  geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('') + ylab('') +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), axis.text.x = element_text(angle = 0), legend.position = "top") + 
  ylab("Impact on Recall of Seeing the Ads") +
  coord_flip()
ggsave(paste0("application/graphs/compareBounds_KB_2021_recall.pdf"), coefs_p, device = "pdf", height = 9, width = 9) 

