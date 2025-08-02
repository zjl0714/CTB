aggBounds_binary_oneshot <- function(dat_bo, seed = NULL, Y, X, S, D, W, Pscore,
                                     splitting = F, cv_fold = cv_fold, trim_l, trim_u){

  scores_l_k <- scores_u_k <- qar_k <- c()
  score1s_l_k <- score1s_u_k <- score0s_k <- score0s_l_k <- score0s_u_k <- score1s_k <- c()
  directions_k <- cor1s_k <- cor0s_k <- qar_help_k <- qar_hurt_k <- c()

  for (k in 1:cv_fold){

    dat.s1 <- dat_bo[dat_bo$split_index != k, ]
    dat.s2 <- dat_bo[dat_bo$split_index == k, ]

    X_eval_k <- dat.s2[, X]
    Y_eval_k <- dat.s2[, Y]
    D_eval_k <- dat.s2[, D]
    S_eval_k <- dat.s2[, S]
    W_eval_k <- dat.s2[, W]
    Ps_eval_k <- dat.s2[, Pscore]

    dat.s1.0 <- dat.s1[dat.s1[, D] == 0, ]
    dat.s1.1 <- dat.s1[dat.s1[, D] == 1, ]

    if (is.null(seed)){
      pf0 <- probability_forest(as.matrix(dat.s1.0[, X]), as.factor(dat.s1.0[, S]), sample.weights = dat.s1.0[, W])
      pf1 <- probability_forest(as.matrix(dat.s1.1[, X]), as.factor(dat.s1.1[, S]), sample.weights = dat.s1.1[, W])
    }else{
      pf0 <- probability_forest(as.matrix(dat.s1.0[, X]), as.factor(dat.s1.0[, S]), sample.weights = dat.s1.0[, W], seed = seed)
      pf1 <- probability_forest(as.matrix(dat.s1.1[, X]), as.factor(dat.s1.1[, S]), sample.weights = dat.s1.1[, W], seed = seed)
    }

    pf0_pred <- predict(pf0, newdata = (as.matrix(X_eval_k)))
    pf1_pred <- predict(pf1, newdata = (as.matrix(X_eval_k)))

    one_q0 <- pf0_pred$predictions[, 2]
    one_q1 <- pf1_pred$predictions[, 2]

    directions <- as.numeric(one_q0 < one_q1)

    one_q_l <- directions * (one_q0 / one_q1) + (1 - directions) * (one_q1 / one_q0)

    one_q_l[one_q_l < trim_l] <- trim_l
    one_q_l[one_q_l > trim_u] <- trim_u

    one_q_u <- 1 - one_q_l

    one_q_l_order <- order(one_q_l)
    one_q_u_order <- order(one_q_u)

    one_q_l_ordered <- one_q_l[one_q_l_order]
    one_q_u_ordered <- one_q_u[one_q_u_order]

    if (is.null(seed)){
      xi_model_help <- probability_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), as.factor(dat.s1.1[dat.s1.1[, S] == 1, Y]))
      xi_model_hurt <- probability_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), as.factor(dat.s1.0[dat.s1.0[, S] == 1, Y]))
    }else{
      xi_model_help <- probability_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), as.factor(dat.s1.1[dat.s1.1[, S] == 1, Y]), seed = seed)
      xi_model_hurt <- probability_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), as.factor(dat.s1.0[dat.s1.0[, S] == 1, Y]), seed = seed)
    }

    xi_pred_help <- predict(xi_model_help, newdata = (as.matrix(X_eval_k)))
    xi_pred_hurt <- predict(xi_model_hurt, newdata = (as.matrix(X_eval_k)))

    xi_help <- xi_pred_help$predictions[, 1]
    xi_hurt <- xi_pred_hurt$predictions[, 1]

    Y_eval_k[is.na(Y_eval_k)] <- 0

    score1s_l_help <- (one_q_l - xi_help) * as.numeric(one_q_l >= xi_help) * one_q1
    score1s_u_help <- (one_q_l + (1 - xi_help - one_q_l) * as.numeric(one_q_l >= 1 - xi_help)) * one_q1

    score0s_help <- (Y_eval_k * S_eval_k) * ((1 - D_eval_k) / (1 - Ps_eval_k))

    score0s_l_hurt <- (one_q_l - xi_hurt) * as.numeric(one_q_l >= xi_hurt) * one_q0
    score0s_u_hurt <- (one_q_l + (1 - xi_hurt - one_q_l) * as.numeric(one_q_l >= 1 - xi_hurt)) * one_q0

    score1s_hurt <- (Y_eval_k * S_eval_k) * (D_eval_k / Ps_eval_k)

    qar_scores_help <- unlist(dat.s2[, S] * (1 - dat.s2[, D]) / (1 - Ps_eval_k))

    qar_scores_hurt <- unlist(dat.s2[, S] * dat.s2[, D] / Ps_eval_k)

    qar_total <- directions * qar_scores_help + (1 - directions) * qar_scores_hurt

    scores_l_total <- directions * (score1s_l_help - score0s_help) + (1 - directions) * (score1s_hurt - score0s_u_hurt)
    scores_u_total <- directions * (score1s_u_help - score0s_help) + (1 - directions) * (score1s_hurt - score0s_l_hurt)

    qar_k <- c(qar_k, qar_total)

    scores_l_k <- c(scores_l_k, scores_l_total)
    scores_u_k <- c(scores_u_k, scores_u_total)
    cat(k, " ")
  }

  # qar_k <- directions_k * qar_help_k + (1 - directions_k) * qar_hurt_k
  # scores_l_k <- directions_k * (score1s_l_k - score0s_k) + (1 - directions_k) * (score1s_k - score0s_u_k)
  # scores_u_k <- directions_k * (score1s_u_k - score0s_k) + (1 - directions_k) * (score1s_k - score0s_l_k)

  N <- sum(dat_bo[, W])
  tau_l_est_avg <- mean(scores_l_k) / mean(qar_k)
  tau_u_est_avg <- mean(scores_u_k) / mean(qar_k)

  V_l <- matrix(NA, 2, 2)
  V_l[1, 1] <- var(scores_l_k)
  V_l[1, 2] <- V_l[2, 1] <- cov(scores_l_k, qar_k)
  V_l[2, 2] <- var(qar_k)

  V_u <- matrix(NA, 2, 2)
  V_u[1, 1] <- var(scores_u_k)
  V_u[1, 2] <- V_u[2, 1] <- cov(scores_u_k, qar_k)
  V_u[2, 2] <- var(qar_k)

  deriv_l <- c(1/mean(qar_k), -tau_l_est_avg/mean(qar_k))
  deriv_u <- c(1/mean(qar_k), -tau_u_est_avg/mean(qar_k))

  se_tau_l_est_avg <- sqrt(t(deriv_l) %*% V_l %*% deriv_l / N)
  se_tau_u_est_avg <- sqrt(t(deriv_u) %*% V_u %*% deriv_u / N)

  result <- c(tau_l_est_avg, tau_u_est_avg, se_tau_l_est_avg, se_tau_u_est_avg)
  return(result)
}

aggBounds_binary <- function(dat_final, seed = NULL, Y, X, S, D, W, Pscore,
                             splitting = F, cv_fold = cv_fold, trim_l, trim_u){

  ests <- aggBounds_binary_oneshot(dat_bo = dat_final, seed = seed, Y, X, S, D, W, Pscore,
                                   splitting = F, cv_fold = cv_fold, trim_l, trim_u)

  tau_l_est_avg <- ests[1]
  tau_u_est_avg <- ests[2]
  se_tau_l_est_avg <- ests[3]
  se_tau_u_est_avg <- ests[4]

  # ndat <- nrow(dat_final)
  # ests_boot <- matrix(NA, Nboots, 2)
  # for (b in 1:Nboots){
  #   dat_boot <- dat_final[sample(ndat, ndat, replace = 1), ]
  #   ests_boot[b, ] <- aggBounds_binary_oneshot(dat_boot, seed = NULL, Y, X, S, D, W, Pscore,
  #                                              regression.splitting = F, cv_fold = cv_fold, trim_l, trim_u)
  # }
  #
  # se_tau_l_est_avg <- sqrt(var(ests_boot[, 1]))
  # se_tau_u_est_avg <- sqrt(var(ests_boot[, 2]))

  result <- list("tau_l_est_avg" = tau_l_est_avg,
                 "tau_u_est_avg" = tau_u_est_avg,
                 "se_tau_l_est_avg" = se_tau_l_est_avg,
                 "se_tau_u_est_avg" = se_tau_u_est_avg)
  return(result)
}


LBounds_binary_oneshot <- function(dat_bo, Y, X, S, D, W, direction){

  q0_lee <- sum(dat_bo[, S] * (1-dat_bo[, D]) * dat_bo[, W]) / sum((1-dat_bo[, D]) * dat_bo[, W])
  q1_lee <- sum(dat_bo[, S] * dat_bo[, D] * dat_bo[, W]) / sum(dat_bo[, D] * dat_bo[, W])
  direction <- as.numeric(q0_lee < q1_lee)

  if (direction){
    q_lee <- q0_lee / q1_lee
    xi_lee <- sum(dat_bo[, S] * dat_bo[, D] * dat_bo[, W] * (1 - dat_bo[, Y]), na.rm = 1) / sum(dat_bo[, S] * dat_bo[, D] * dat_bo[, W])
    theta1_l_lee <- (q_lee - xi_lee) * (q_lee >= xi_lee) / q_lee
    theta1_u_lee <- (q_lee + (1 - q_lee - xi_lee) * (1 - q_lee <= xi_lee)) / q_lee
    theta0_lee <- sum((1-dat_bo[, D]) * dat_bo[, S] * dat_bo[, W] * dat_bo[, Y], na.rm = 1) / sum((1-dat_bo[, D]) * dat_bo[, S] * dat_bo[, W], na.rm = 1)
    tau_l_lee <- theta1_l_lee - theta0_lee
    tau_u_lee <- theta1_u_lee - theta0_lee
  }else{
    q_lee <- q1_lee / q0_lee
    xi_lee <- sum(dat_bo[, S] * (1 - dat_bo[, D]) * dat_bo[, W] * (1 - dat_bo[, Y]), na.rm = 1) / sum(dat_bo[, S] * (1 - dat_bo[, D]) * dat_bo[, W])
    theta0_l_lee <- (q_lee - xi_lee) * (q_lee >= xi_lee) / q_lee
    theta0_u_lee <- (q_lee + (1 - q_lee - xi_lee) * (1 - q_lee <= xi_lee)) / q_lee
    theta1_lee <- sum(dat_bo[, D] * dat_bo[, S] * dat_bo[, W] * dat_bo[, Y], na.rm = 1) / sum(dat_bo[, D] * dat_bo[, S] * dat_bo[, W], na.rm = 1)
    tau_l_lee <- theta1_lee - theta0_u_lee
    tau_u_lee <- theta1_lee - theta0_l_lee
  }
  return(c(tau_l_lee, tau_u_lee))
}


LBounds_binary <- function(dat_final, Y, X, S, D, W, Pscore,
                           cond.mono = FALSE, Nboots = 1000){

  if (cond.mono){
    if (is.null(seed)){
      pf0 <- probability_forest(X = as.matrix(dat_final[dat_final[, D] == 0, X]), Y = as.factor(unlist(dat_final[dat_final[, D] == 0, S])), sample.weights = dat_final[dat_final[, D] == 0, W])
      pf1 <- probability_forest(X = as.matrix(dat_final[dat_final[, D] == 1, X]), Y = as.factor(unlist(dat_final[dat_final[, D] == 1, S])), sample.weights = dat_final[dat_final[, D] == 1, W])
    }else{
      pf0 <- probability_forest(X = as.matrix(dat_final[dat_final[, D] == 0, X]), Y = as.factor(unlist(dat_final[dat_final[, D] == 0, S])), sample.weights = dat_final[dat_final[, D] == 0, W], seed = seed)
      pf1 <- probability_forest(X = as.matrix(dat_final[dat_final[, D] == 1, X]), Y = as.factor(unlist(dat_final[dat_final[, D] == 1, S])), sample.weights = dat_final[dat_final[, D] == 1, W], seed = seed)
    }

    pf0_pred <- predict(pf0, newdata = (as.matrix(dat_final[, X])))
    pf1_pred <- predict(pf1, newdata = (as.matrix(dat_final[, X])))

    one_q0 <- pf0_pred$predictions[, 2]
    one_q1 <- pf1_pred$predictions[, 2]

    directions <- as.numeric(one_q0 < one_q1)

    dat_help <- dat_final[directions == 1, ]
    dat_hurt <- dat_final[directions == 0, ]

    ndat_help <- nrow(dat_help)
    ndat_hurt <- nrow(dat_hurt)

    # hurt lee bounds
    ests_lee_hurt <- LBounds_binary_oneshot(dat_hurt, Y, X, S, D, W, direction = 0)

    tau_u_lee_hurt <- ests_lee_hurt[2]
    tau_l_lee_hurt <- ests_lee_hurt[1]

    ests_lee_boot_hurt <- matrix(NA, Nboots, 2)
    for (b in 1:Nboots){
      dat_hurt_boot <- dat_hurt[sample(ndat_hurt, ndat_hurt, replace = 1), ]
      ests_lee_boot_hurt[b, ] <- LBounds_binary_oneshot(dat_hurt_boot, Y, X, S, D, W, direction = 0)
      if (b %% 100 == 0){
        cat(b, "\n")
      }
    }

    se_tau_l_lee_hurt <- sqrt(var(ests_lee_boot_hurt[, 1]))
    se_tau_u_lee_hurt <- sqrt(var(ests_lee_boot_hurt[, 2]))

    # help lee bounds
    ests_lee_help <- LBounds_binary_oneshot(dat_help, Y, X, S, D, W, direction = 1)

    tau_u_lee_help <- ests_lee_help[2]
    tau_l_lee_help <- ests_lee_help[1]

    ests_lee_boot_help <- matrix(NA, Nboots, 2)
    for (b in 1:Nboots){
      dat_help_boot <- dat_help[sample(ndat_help, ndat_help, replace = 1), ]
      ests_lee_boot_help[b, ] <- LBounds_binary_oneshot(dat_help_boot, Y, X, S, D, W, direction = 1)
      if (b %% 100 == 0){
        cat(b, "\n")
      }
    }

    se_tau_l_lee_help <- sqrt(var(ests_lee_boot_help[, 1]))
    se_tau_u_lee_help <- sqrt(var(ests_lee_boot_help[, 2]))

    pos_ratio <- mean(directions)

    tau_u_lee <- pos_ratio * tau_u_lee_help + (1-pos_ratio) * tau_u_lee_hurt
    tau_l_lee <- pos_ratio * tau_l_lee_help + (1-pos_ratio) * tau_l_lee_hurt

    se_tau_l_lee <- sqrt(pos_ratio^2 * se_tau_l_lee_help^2 + (1-pos_ratio)^2 * se_tau_l_lee_hurt^2)
    se_tau_u_lee <- sqrt(pos_ratio^2 * se_tau_u_lee_help^2 + (1-pos_ratio)^2 * se_tau_u_lee_hurt^2)

  }else{
    ndat <- nrow(dat_final)
    q0_lee <- sum(dat_final[, S] * (1-dat_final[, D]) * dat_final[, W]) / sum((1-dat_final[, D]) * dat_final[, W])
    q1_lee <- sum(dat_final[, S] * dat_final[, D] * dat_final[, W]) / sum(dat_final[, D] * dat_final[, W])

    direction <- as.numeric(q0_lee < q1_lee)
    if (direction){
      cat("S(1) > S(0)", "\n")

      ests_lee <- LBounds_binary_oneshot(dat_final, Y, X, S, D, W, direction = 1)

      tau_u_lee <- ests_lee[2]
      tau_l_lee <- ests_lee[1]

      ests_lee_boot <- matrix(NA, Nboots, 2)
      for (b in 1:Nboots){
        dat_final_boot <- dat_final[sample(ndat, ndat, replace = 1), ]
        ests_lee_boot[b, ] <- LBounds_binary_oneshot(dat_final_boot, Y, X, S, D, W, direction = 1)
      }

      se_tau_l_lee <- sqrt(var(ests_lee_boot[, 1]))
      se_tau_u_lee <- sqrt(var(ests_lee_boot[, 2]))
    }else{
      cat("S(0) > S(1)", "\n")
      ests_lee <- LBounds_binary_oneshot(dat_final, Y, X, S, D, W, direction = 0)

      tau_u_lee <- ests_lee[2]
      tau_l_lee <- ests_lee[1]

      ests_lee_boot <- matrix(NA, Nboots, 2)
      for (b in 1:Nboots){
        dat_final_boot <- dat_final[sample(ndat, ndat, replace = 1), ]
        ests_lee_boot[b, ] <- LBounds_binary_oneshot(dat_final_boot, Y, X, S, D, W, direction = 0)
      }

      se_tau_l_lee <- sqrt(var(ests_lee_boot[, 1]))
      se_tau_u_lee <- sqrt(var(ests_lee_boot[, 2]))
    }
  }
  result <- list("tau_l_est_lee" = tau_l_lee,
                 "tau_u_est_lee" = tau_u_lee,
                 "se_tau_l_est_lee" = se_tau_l_lee,
                 "se_tau_u_est_lee" = se_tau_u_lee)
  return(result)
}
