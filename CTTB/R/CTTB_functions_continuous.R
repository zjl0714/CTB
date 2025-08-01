# functions for continuous outcomes

aggBounds_cont <- function(dat_final, seed = NULL, Y, X, S, D, W, Pscore,
                           splitting = F, cv_fold = cv_fold, trim_l, trim_u){

  scores_l_k <- scores_u_k <- qar_k <- Ps_k <- c()
  score1s_l_k <- score1s_u_k <- score0s_k <- score0s_l_k <- score0s_u_k <- score1s_k <- c()
  directions_k <- cor1s_k <- cor0s_k <- qar_help_k <- qar_hurt_k <- c()

  for (k in 1:cv_fold){

    dat.s1 <- dat_final[dat_final$split_index != k, ]
    dat.s2 <- dat_final[dat_final$split_index == k, ]

    X_eval_k <- dat.s2[, X]
    Y_eval_k <- dat.s2[, Y]
    D_eval_k <- dat.s2[, D]
    S_eval_k <- dat.s2[, S]
    W_eval_k <- dat.s2[, W]
    if (is.null(Pscore)){
      s_index <- sample(nrow(dat.s1), nrow(dat.s1)/2)
      dat.s3 <- dat.s1[s_index, ]
      dat.s1 <- dat.s1[-s_index, ]
      if (is.null(seed)){
        ps <- probability_forest(as.matrix(dat.s3[, X]), as.factor(dat.s3[, D]), sample.weights = dat.s3[, W])
      }else{
        ps <- probability_forest(as.matrix(dat.s3[, X]), as.factor(dat.s3[, D]), sample.weights = dat.s3[, W], seed = seed)
      }
      Ps_pred <- predict(ps, newdata = (as.matrix(X_eval_k)))
      Ps_eval_k <- Ps_pred$predictions[, 2]
      Ps_eval_k[Ps_eval_k <= trim_l] <- trim_l
      Ps_eval_k[Ps_eval_k >= trim_u] <- trim_u
    }else{
      Ps_eval_k <- dat.s2[, Pscore]
    }

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

    directions <- as.numeric(one_q0 / one_q1 < (1 + trim_l))

    one_q_l <- directions * (one_q0 / one_q1) + (1 - directions) * (one_q1 / one_q0)

    one_q_l[one_q_l < trim_l] <- trim_l
    one_q_l[one_q_l > trim_u] <- trim_u

    one_q_u <- 1 - one_q_l

    one_q_l_order <- order(one_q_l)
    one_q_u_order <- order(one_q_u)

    one_q_l_ordered <- one_q_l[one_q_l_order]
    one_q_u_ordered <- one_q_u[one_q_u_order]

    if (is.null(seed)){
      q1.forest.lower <- quantile_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), dat.s1.1[dat.s1.1[, S] == 1, Y], quantiles = one_q_l_ordered, regression.splitting = splitting)
      q1.forest.upper <- quantile_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), dat.s1.1[dat.s1.1[, S] == 1, Y], quantiles = one_q_u_ordered, regression.splitting = splitting)

      q0.forest.lower <- quantile_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), dat.s1.0[dat.s1.0[, S] == 1, Y], quantiles = one_q_l_ordered, regression.splitting = splitting)
      q0.forest.upper <- quantile_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), dat.s1.0[dat.s1.0[, S] == 1, Y], quantiles = one_q_u_ordered, regression.splitting = splitting)
    }else{
      q1.forest.lower <- quantile_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), dat.s1.1[dat.s1.1[, S] == 1, Y], seed = seed, quantiles = one_q_l_ordered, regression.splitting = splitting)
      q1.forest.upper <- quantile_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), dat.s1.1[dat.s1.1[, S] == 1, Y], seed = seed, quantiles = one_q_u_ordered, regression.splitting = splitting)

      q0.forest.lower <- quantile_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), dat.s1.0[dat.s1.0[, S] == 1, Y], seed = seed, quantiles = one_q_l_ordered, regression.splitting = splitting)
      q0.forest.upper <- quantile_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), dat.s1.0[dat.s1.0[, S] == 1, Y], seed = seed, quantiles = one_q_u_ordered, regression.splitting = splitting)
    }

    q1.forest.lower.pred <- predict(q1.forest.lower, newdata = (as.matrix(X_eval_k)), quantiles = one_q_l_ordered)
    q1.forest.upper.pred <- predict(q1.forest.upper, newdata = (as.matrix(X_eval_k)), quantiles = one_q_u_ordered)

    q0.forest.lower.pred <- predict(q0.forest.lower, newdata = (as.matrix(X_eval_k)), quantiles = one_q_l_ordered)
    q0.forest.upper.pred <- predict(q0.forest.upper, newdata = (as.matrix(X_eval_k)), quantiles = one_q_u_ordered)

    yq_l_help <- sapply(1:length(one_q_l), function(tt){
      return(q1.forest.lower.pred$predictions[tt, which(one_q_l_order==tt)])
    })
    yq_u_help <- sapply(1:length(one_q_u), function(tt){
      return(q1.forest.upper.pred$predictions[tt, which(one_q_u_order==tt)])
    })

    yq_l_hurt <- sapply(1:length(one_q_l), function(tt){
      return(q0.forest.lower.pred$predictions[tt, which(one_q_l_order==tt)])
    })
    yq_u_hurt <- sapply(1:length(one_q_u), function(tt){
      return(q0.forest.upper.pred$predictions[tt, which(one_q_u_order==tt)])
    })

    yq_l <- directions*yq_l_help + (1-directions)*yq_l_hurt
    yq_u <- directions*yq_u_help + (1-directions)*yq_u_hurt


    Y_eval_k[is.na(Y_eval_k)] <- 0

    score1s_l <- ((Y_eval_k - yq_l) * as.numeric(Y_eval_k <= yq_l) * S_eval_k * (D_eval_k / Ps_eval_k) +
      S_eval_k * yq_l * (1 - D_eval_k) / (1 - Ps_eval_k)) * W_eval_k
    score1s_u <- ((Y_eval_k - yq_u) * as.numeric(Y_eval_k >= yq_u) * S_eval_k * (D_eval_k / Ps_eval_k) +
      S_eval_k * yq_u * (1 - D_eval_k) / (1 - Ps_eval_k)) * W_eval_k

    score0s <- (Y_eval_k * S_eval_k * (1 - D_eval_k) / (1 - Ps_eval_k)) * W_eval_k

    score0s_l <- ((Y_eval_k - yq_l) * as.numeric(Y_eval_k <= yq_l) * S_eval_k * (1 - D_eval_k) / (1 - Ps_eval_k) +
      S_eval_k * yq_l * D_eval_k / Ps_eval_k) * W_eval_k
    score0s_u <- ((Y_eval_k - yq_u) * as.numeric(Y_eval_k >= yq_u) * S_eval_k * (1 - D_eval_k) / (1 - Ps_eval_k) +
      S_eval_k * yq_u * D_eval_k / Ps_eval_k) * W_eval_k

    score1s <- (Y_eval_k * S_eval_k * D_eval_k / Ps_eval_k) * W_eval_k

    qar_scores_help <- S_eval_k * (1 - D_eval_k) / (1 - Ps_eval_k) * W_eval_k
    qar_scores_hurt <- S_eval_k * D_eval_k / Ps_eval_k * W_eval_k

    if (is.null(Pscore)){
      # score1s_l_help <- score1s_l_help + ((yq_l * one_q0 / (Ps_eval_k * (1 - Ps_eval_k))) - (sum(score1s_l_help) / sum(W_eval_k)) / Ps_eval_k) * (D_eval_k - Ps_eval_k) * W_eval_k
      # score1s_u_help <- score1s_u_help + ((yq_u * one_q0 / (Ps_eval_k * (1 - Ps_eval_k))) - (sum(score1s_u_help) / sum(W_eval_k)) / Ps_eval_k) * (D_eval_k - Ps_eval_k) * W_eval_k
      # score0s_help <- score0s_help + (sum(score0s_help) / sum(W_eval_k)) / (1 - Ps_eval_k) * (D_eval_k - Ps_eval_k) * W_eval_k
      # score0s_l_hurt <- score0s_l_hurt - ((yq_l * one_q0 / (Ps_eval_k * (1 - Ps_eval_k))) - (sum(score0s_l_hurt) / sum(W_eval_k)) / (1 - Ps_eval_k)) * (D_eval_k - Ps_eval_k) * W_eval_k
      # score0s_u_hurt <- score0s_u_hurt - ((yq_u * one_q0 / (Ps_eval_k * (1 - Ps_eval_k))) - (sum(score0s_u_hurt) / sum(W_eval_k)) / (1 - Ps_eval_k)) * (D_eval_k - Ps_eval_k) * W_eval_k
      # score1s_hurt <- score1s_hurt - (sum(score1s_hurt) / sum(W_eval_k)) / Ps_eval_k * (D_eval_k - Ps_eval_k) * W_eval_k
      score1s_l <- score1s_l + (yq_l * one_q0 * (D_eval_k / Ps_eval_k - (1 - D_eval_k) / (1 - Ps_eval_k)) * W_eval_k)
      score1s_u <- score1s_u + (yq_u * one_q0 * (D_eval_k / Ps_eval_k - (1 - D_eval_k) / (1 - Ps_eval_k)) * W_eval_k)
      score0s_l <- score0s_l - (yq_l * one_q1 * (D_eval_k / Ps_eval_k - (1 - D_eval_k) / (1 - Ps_eval_k)) * W_eval_k)
      score0s_u <- score0s_u - (yq_u * one_q1 * (D_eval_k / Ps_eval_k - (1 - D_eval_k) / (1 - Ps_eval_k)) * W_eval_k)
      cor1s_k <- c(cor1s_k, D_eval_k * W_eval_k / Ps_eval_k)
      cor0s_k <- c(cor0s_k, (1 - D_eval_k) * W_eval_k / (1 - Ps_eval_k))
    }

    score1s_l_k <- c(score1s_l_k, score1s_l)
    score1s_u_k <- c(score1s_u_k, score1s_u)
    score0s_k <- c(score0s_k, score0s)
    score0s_l_k <- c(score0s_l_k, score0s_l)
    score0s_u_k <- c(score0s_u_k, score0s_u)
    score1s_k <- c(score1s_k, score1s)
    Ps_k <- c(Ps_k, Ps_eval_k)

    qar_help_k <- c(qar_help_k, qar_scores_help)
    qar_hurt_k <- c(qar_hurt_k, qar_scores_hurt)

    directions_k <- c(directions_k, directions)
    # cat(k, " ")
  }
  if (is.null(Pscore)){
    scores_dat <- data.frame("score1s_l_k" = score1s_l_k, "score1s_u_k" = score1s_u_k,
                             "score0s_l_k" = score0s_l_k, "score0s_u_k" = score0s_u_k,
                             "score0s_k" = score0s_k, "score1s_k" = score1s_k,
                             "qar_help_k" = qar_help_k, "qar_hurt_k" = qar_hurt_k,
                             "directions_k" = directions_k, "cor1s_k" = cor1s_k,
                             "cor0s_k" = cor0s_k, "weight" = dat_final[, W])

    N <- sum(scores_dat[, "weight"])
    N_help <- N_hurt <- 0
    tau_l_est_avg_help <- tau_u_est_avg_help <- tau_l_est_avg_hurt <- tau_u_est_avg_hurt <- 0
    var_tau_l_est_avg_help <- var_tau_u_est_avg_help <- var_tau_l_est_avg_hurt <- var_tau_u_est_avg_hurt <- 0
    if (1 %in% directions_k & mean(scores_dat[scores_dat[, "directions_k"] == 1, "cor1s_k"]) > 0 & mean(scores_dat[scores_dat[, "directions_k"] == 1, "cor0s_k"]) > 0 & mean(scores_dat[scores_dat[, "directions_k"] == 1, "qar_help_k"]) > 0){
      scores_dat_help <- scores_dat[scores_dat[, "directions_k"] == 1, ]
      N_help <- sum(scores_dat_help[, "weight"])
      tau_l_est_avg_help <- (mean(scores_dat_help$score1s_l_k) / mean(scores_dat_help$cor1s_k) - mean(scores_dat_help$score0s_k) / mean(scores_dat_help$cor0s_k)) / (mean(scores_dat_help$qar_help_k) / mean(scores_dat_help$cor0s_k))
      tau_u_est_avg_help <- (mean(scores_dat_help$score1s_u_k) / mean(scores_dat_help$cor1s_k) - mean(scores_dat_help$score0s_k) / mean(scores_dat_help$cor0s_k)) / (mean(scores_dat_help$qar_help_k) / mean(scores_dat_help$cor0s_k))

      scores_l_help_k <- scores_dat_help$score1s_l_k / mean(scores_dat_help$cor1s_k) - mean(scores_dat_help$score1s_l_k) * scores_dat_help$cor1s_k / (mean(scores_dat_help$cor1s_k))^2 -
        scores_dat_help$score0s_k / mean(scores_dat_help$cor0s_k) + mean(scores_dat_help$score0s_k) * scores_dat_help$cor0s_k / (mean(scores_dat_help$cor0s_k))^2
      scores_u_help_k <- scores_dat_help$score1s_u_k / mean(scores_dat_help$cor1s_k) - mean(scores_dat_help$score1s_u_k) * scores_dat_help$cor1s_k / (mean(scores_dat_help$cor1s_k))^2 -
        scores_dat_help$score0s_k / mean(scores_dat_help$cor0s_k) + mean(scores_dat_help$score0s_k) * scores_dat_help$cor0s_k / (mean(scores_dat_help$cor0s_k))^2

      qar_help_k <- scores_dat_help$qar_help_k / mean(scores_dat_help$cor0s_k) - mean(scores_dat_help$qar_help_k) * scores_dat_help$cor0s_k / (mean(scores_dat_help$cor0s_k))^2

      V_l_help <- matrix(NA, 2, 2)
      V_l_help[1, 1] <- var(scores_l_help_k)
      V_l_help[1, 2] <- V_l_help[2, 1] <- cov(scores_l_help_k, qar_help_k)
      V_l_help[2, 2] <- var(qar_help_k)

      V_u_help <- matrix(NA, 2, 2)
      V_u_help[1, 1] <- var(scores_u_help_k)
      V_u_help[1, 2] <- V_u_help[2, 1] <- cov(scores_u_help_k, qar_help_k)
      V_u_help[2, 2] <- var(qar_help_k)

      deriv_l_help <- c(mean(scores_dat_help$cor0s_k)/mean(scores_dat_help$qar_help_k), -tau_l_est_avg_help * mean(scores_dat_help$cor0s_k) / mean(scores_dat_help$qar_help_k))
      deriv_u_help <- c(mean(scores_dat_help$cor0s_k)/mean(scores_dat_help$qar_help_k), -tau_u_est_avg_help * mean(scores_dat_help$cor0s_k) / mean(scores_dat_help$qar_help_k))

      var_tau_l_est_avg_help <- t(deriv_l_help) %*% V_l_help %*% deriv_l_help
      var_tau_u_est_avg_help <- t(deriv_u_help) %*% V_u_help %*% deriv_u_help
    }

    if (0 %in% directions_k & mean(scores_dat[scores_dat[, "directions_k"] == 0, "cor1s_k"]) > 0 & mean(scores_dat[scores_dat[, "directions_k"] == 0, "cor0s_k"]) > 0 & mean(scores_dat[scores_dat[, "directions_k"] == 0, "qar_hurt_k"]) > 0){
      scores_dat_hurt <- scores_dat[scores_dat[, "directions_k"] == 0, ]
      N_hurt <- sum(scores_dat_hurt[, "weight"])
      tau_l_est_avg_hurt <- (mean(scores_dat_hurt$score1s_k) / mean(scores_dat_hurt$cor1s_k) - mean(scores_dat_hurt$score0s_u_k) / mean(scores_dat_hurt$cor0s_k)) / (mean(scores_dat_hurt$qar_hurt_k) / mean(scores_dat_hurt$cor1s_k))
      tau_u_est_avg_hurt <- (mean(scores_dat_hurt$score1s_k) / mean(scores_dat_hurt$cor1s_k) - mean(scores_dat_hurt$score0s_l_k) / mean(scores_dat_hurt$cor0s_k)) / (mean(scores_dat_hurt$qar_hurt_k) / mean(scores_dat_hurt$cor1s_k))

      scores_l_hurt_k <- scores_dat_hurt$score1s_k / mean(scores_dat_hurt$cor1s_k) - mean(scores_dat_hurt$score1s_k) * scores_dat_hurt$cor1s_k / (mean(scores_dat_hurt$cor1s_k))^2 -
        scores_dat_hurt$score0s_u_k / mean(scores_dat_hurt$cor0s_k) + mean(scores_dat_hurt$score0s_u_k) * scores_dat_hurt$cor0s_k / (mean(scores_dat_hurt$cor0s_k))^2
      scores_u_hurt_k <- scores_dat_hurt$score1s_k / mean(scores_dat_hurt$cor1s_k) - mean(scores_dat_hurt$score1s_k) * scores_dat_hurt$cor1s_k / (mean(scores_dat_hurt$cor1s_k))^2 -
        scores_dat_hurt$score0s_l_k / mean(scores_dat_hurt$cor0s_k) + mean(scores_dat_hurt$score0s_l_k) * scores_dat_hurt$cor0s_k / (mean(scores_dat_hurt$cor0s_k))^2

      qar_hurt_k <- scores_dat_hurt$qar_hurt_k / mean(scores_dat_hurt$cor1s_k) - mean(scores_dat_hurt$qar_hurt_k) * scores_dat_hurt$cor1s_k / (mean(scores_dat_hurt$cor1s_k))^2

      V_l_hurt <- matrix(NA, 2, 2)
      V_l_hurt[1, 1] <- var(scores_l_hurt_k)
      V_l_hurt[1, 2] <- V_l_hurt[2, 1] <- cov(scores_l_hurt_k, qar_hurt_k)
      V_l_hurt[2, 2] <- var(qar_hurt_k)

      V_u_hurt <- matrix(NA, 2, 2)
      V_u_hurt[1, 1] <- var(scores_u_hurt_k)
      V_u_hurt[1, 2] <- V_u_hurt[2, 1] <- cov(scores_u_hurt_k, qar_hurt_k)
      V_u_hurt[2, 2] <- var(qar_hurt_k)

      deriv_l_hurt <- c(mean(scores_dat_hurt$cor1s_k) / mean(scores_dat_hurt$qar_hurt_k), -tau_l_est_avg_hurt * mean(scores_dat_hurt$cor1s_k) / mean(scores_dat_hurt$qar_hurt_k))
      deriv_u_hurt <- c(mean(scores_dat_hurt$cor1s_k) / mean(scores_dat_hurt$qar_hurt_k), -tau_u_est_avg_hurt * mean(scores_dat_hurt$cor1s_k) / mean(scores_dat_hurt$qar_hurt_k))

      var_tau_l_est_avg_hurt <- t(deriv_l_hurt) %*% V_l_hurt %*% deriv_l_hurt
      var_tau_u_est_avg_hurt <- t(deriv_u_hurt) %*% V_u_hurt %*% deriv_u_hurt
    }

    tau_l_est_avg <- (N_help * tau_l_est_avg_help + N_hurt * tau_l_est_avg_hurt) / N
    tau_u_est_avg <- (N_help * tau_u_est_avg_help + N_hurt * tau_u_est_avg_hurt) / N

    se_tau_l_est_avg <- sqrt((var_tau_l_est_avg_help + var_tau_l_est_avg_hurt) / N)
    se_tau_u_est_avg <- sqrt((var_tau_u_est_avg_help + var_tau_u_est_avg_hurt) / N)
  }else{
    qar_k <- directions_k * qar_help_k + (1 - directions_k) * qar_hurt_k
    scores_l_k <- directions_k * (score1s_l_k - score0s_k) + (1 - directions_k) * (score1s_k - score0s_u_k)
    scores_u_k <- directions_k * (score1s_u_k - score0s_k) + (1 - directions_k) * (score1s_k - score0s_l_k)

    N <- sum(dat_final[, W])
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
  }


  result <- list("scores_l" = scores_l_k, "scores_u" = scores_u_k,
                 "trimming_prob" = qar_k, "tau_l_est_avg" = tau_l_est_avg,
                 "tau_u_est_avg" = tau_u_est_avg, "se_tau_l_est_avg" = se_tau_l_est_avg,
                 "se_tau_u_est_avg" = se_tau_u_est_avg, "Ps_k" = Ps_k)
  return(result)
}


condBounds_cont <- function(dat_final, seed=NULL, Y, X, S, D, W, Pscore,
                            X_moderator, splitting = F,
                            cv_fold = cv_fold){

  tau_m <- se_tau_m <- tau_l_m <- se_tau_l_m <- tau_u_m <- se_tau_u_m <- matrix(NA, nrow(X_moderator), cv_fold)

  for (k in 1:cv_fold){
    dat.s1 <- dat_final[dat_final$split_index != k, ]
    dat.s2 <- dat_final[dat_final$split_index == k, ]

    X_eval_k <- dat.s2[, X]
    Y_eval_k <- dat.s2[, Y]
    D_eval_k <- dat.s2[, D]
    S_eval_k <- dat.s2[, S]
    W_eval_k <- dat.s2[, W]
    if (is.null(Pscore)){
      if (is.null(seed)){
        ps <- probability_forest(as.matrix(dat.s1[, X]), as.factor(dat.s1[, D]), sample.weights = dat.s1[, W])
      }else{
        ps <- probability_forest(as.matrix(dat.s1[, X]), as.factor(dat.s1[, D]), sample.weights = dat.s1[, W], seed = seed)
      }
      Ps_pred <- predict(ps, newdata = (as.matrix(X_eval_k)))
      Ps_eval_k <- Ps_pred$predictions[, 2]
    }else{
      Ps_eval_k <- dat.s2[, Pscore]
    }

    dat.s1.0 <- dat.s1[dat.s1[, D] == 0, ]
    dat.s1.1 <- dat.s1[dat.s1[, D] == 1, ]

    if (is.null(seed)){
      pf0 <- probability_forest(as.matrix(dat.s1.0[, X]), as.factor(dat.s1.0[, S]))
      pf1 <- probability_forest(as.matrix(dat.s1.1[, X]), as.factor(dat.s1.1[, S]))
    }else{
      pf0 <- probability_forest(as.matrix(dat.s1.0[, X]), as.factor(dat.s1.0[, S]), seed = seed)
      pf1 <- probability_forest(as.matrix(dat.s1.1[, X]), as.factor(dat.s1.1[, S]), seed = seed)
    }

    pf0_pred_m <- predict(pf0, newdata = (as.matrix(X_moderator)), estimate.variance = 1)
    pf1_pred_m <- predict(pf1, newdata = (as.matrix(X_moderator)), estimate.variance = 1)

    one_q0_m <- pf0_pred_m$predictions[, 2]
    one_q1_m <- pf1_pred_m$predictions[, 2]

    se_one_q0_m <- sqrt(pf0_pred_m$variance.estimates[, 2])
    se_one_q1_m <- sqrt(pf1_pred_m$variance.estimates[, 2])

    directions_m <- as.numeric(one_q0_m < one_q1_m)
    one_q_l_m <- directions_m*(one_q0_m / one_q1_m) + (1-directions_m)*(one_q1_m / one_q0_m)
    se_one_q_l_m <- directions_m*(se_one_q0_m / one_q1_m + se_one_q1_m * (one_q0_m / one_q1_m^2)) +
      (1-directions_m)*(se_one_q1_m / one_q0_m + se_one_q0_m * (one_q1_m / one_q0_m^2))

    one_q_u_m <- 1 - one_q_l_m

    one_q_l_order_m <- order(one_q_l_m)
    one_q_u_order_m <- order(one_q_u_m)

    one_q_l_ordered_m <- one_q_l_m[one_q_l_order_m]
    one_q_u_ordered_m <- one_q_u_m[one_q_u_order_m]

    if (is.null(seed)){
      q1.forest.lower.m <- quantile_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), dat.s1.1[dat.s1.1[, S] == 1, Y], quantiles = one_q_l_ordered_m, regression.splitting = 0)
      q1.forest.upper.m <- quantile_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), dat.s1.1[dat.s1.1[, S] == 1, Y], quantiles = one_q_u_ordered_m, regression.splitting = 0)
    }else{
      q1.forest.lower.m <- quantile_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), dat.s1.1[dat.s1.1[, S] == 1, Y], seed = seed, quantiles = one_q_l_ordered_m, regression.splitting = 0)
      q1.forest.upper.m <- quantile_forest(as.matrix(dat.s1.1[dat.s1.1[, S] == 1, X]), dat.s1.1[dat.s1.1[, S] == 1, Y], seed = seed, quantiles = one_q_u_ordered_m, regression.splitting = 0)
    }
    q1.forest.lower.pred.m <- predict(q1.forest.lower.m, newdata = (as.matrix(X_moderator)), quantiles = one_q_l_ordered_m)
    q1.forest.upper.pred.m <- predict(q1.forest.upper.m, newdata = (as.matrix(X_moderator)), quantiles = one_q_u_ordered_m)


    yq_l_help_m <- sapply(1:length(one_q_l_m), function(tt){
      return(q1.forest.lower.pred.m$predictions[tt, which(one_q_l_order_m==tt)])
    })
    yq_u_help_m <- sapply(1:length(one_q_l_m), function(tt){
      return(q1.forest.upper.pred.m$predictions[tt, which(one_q_u_order_m==tt)])
    })

    if (is.null(seed)){
      q0.forest.lower.m <- quantile_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), dat.s1.0[dat.s1.0[, S] == 1, Y], quantiles = one_q_l_ordered_m, regression.splitting = 0)
      q0.forest.upper.m <- quantile_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), dat.s1.0[dat.s1.0[, S] == 1, Y], quantiles = one_q_u_ordered_m, regression.splitting = 0)
    }else{
      q0.forest.lower.m <- quantile_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), dat.s1.0[dat.s1.0[, S] == 1, Y], seed = seed, quantiles = one_q_l_ordered_m, regression.splitting = 0)
      q0.forest.upper.m <- quantile_forest(as.matrix(dat.s1.0[dat.s1.0[, S] == 1, X]), dat.s1.0[dat.s1.0[, S] == 1, Y], seed = seed, quantiles = one_q_u_ordered_m, regression.splitting = 0)
    }
    q0.forest.lower.pred.m <- predict(q0.forest.lower.m, newdata = (as.matrix(X_moderator)), quantiles = one_q_l_ordered_m)
    q0.forest.upper.pred.m <- predict(q0.forest.upper.m, newdata = (as.matrix(X_moderator)), quantiles = one_q_u_ordered_m)

    yq_l_hurt_m <- sapply(1:length(one_q_l_m), function(tt){
      return(q0.forest.lower.pred.m$predictions[tt, which(one_q_l_order_m==tt)])
    })
    yq_u_hurt_m <- sapply(1:length(one_q_l_m), function(tt){
      return(q0.forest.upper.pred.m$predictions[tt, which(one_q_u_order_m==tt)])
    })

    yq_l_m <- directions_m*yq_l_help_m + (1-directions_m)*yq_l_hurt_m
    yq_u_m <- directions_m*yq_u_help_m + (1-directions_m)*yq_u_hurt_m

    if (is.null(seed)){
      rf0_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, X]), Y = dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, Y])
    }else{
      rf0_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, X]), Y = dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, Y], seed = seed)
    }
    pred0_m <- predict(rf0_m, newdata = as.matrix(X_moderator), estimate.variance = 1)

    theta0_m <- pred0_m[, 1]
    se_theta0_m <- sqrt(pred0_m[, 2])

    if (is.null(seed)){
      rf1_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, X]), Y = dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, Y])
    }else{
      rf1_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, X]), Y = dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, Y], seed = seed)
    }
    pred1_m <- predict(rf1_m, newdata = as.matrix(X_moderator), estimate.variance = 1)

    theta1_m <- pred1_m[, 1]
    se_theta1_m <- sqrt(pred1_m[, 2])

    if (is.null(seed)){
      cf_m <- causal_forest(X = as.matrix(dat.s2[dat.s2[, S] == 1, X]),
                            Y = dat.s2[dat.s2[, S] == 1, Y],
                            W = dat.s2[dat.s2[, S] == 1, D])
    }else{
      cf_m <- causal_forest(X = as.matrix(dat.s2[dat.s2[, S] == 1, X]),
                            Y = dat.s2[dat.s2[, S] == 1, Y],
                            W = dat.s2[dat.s2[, S] == 1, D], seed = seed)
    }

    pred_m <- predict(cf_m, newdata = as.matrix(X_moderator), estimate.variance = 1)
    tau_m[, k] <- pred_m[, 1]
    se_tau_m[, k] <- pred_m[, 2]

    for (m in 1:nrow(X_moderator)){
      Y1.trim.lower <- (dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, Y] - yq_l_m[m]) * (dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, Y] <= yq_l_m[m])
      Y1.trim.upper <- (dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, Y] - yq_u_m[m]) * (dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, Y] >= yq_u_m[m])
      if (is.null(seed)){
        rf1_l_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, X]), Y = Y1.trim.lower)
        rf1_u_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, X]), Y = Y1.trim.upper)
      }else{
        rf1_l_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, X]), Y = Y1.trim.lower, seed = seed)
        rf1_u_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 1 & dat.s2[, S] == 1, X]), Y = Y1.trim.upper, seed = seed)
      }
      pred1_l_m <- predict(rf1_l_m, newdata = as.matrix(X_moderator), estimate.variance = 1)
      pred1_u_m <- predict(rf1_u_m, newdata = as.matrix(X_moderator), estimate.variance = 1)
      theta1_l_m <- pred1_l_m[m, 1] / one_q_l_m[m] + yq_l_m[m]
      theta1_u_m <- pred1_u_m[m, 1] / one_q_l_m[m] + yq_u_m[m]

      se_theta1_l_m <- sqrt(pred1_l_m[m, 2]) / one_q_l_m[m] + se_one_q_l_m[m] * (pred1_l_m[m, 1] / one_q_l_m[m]^2)
      se_theta1_u_m <- sqrt(pred1_u_m[m, 2]) / one_q_l_m[m] + se_one_q_l_m[m] * (pred1_u_m[m, 1] / one_q_l_m[m]^2)

      Y0.trim.lower <- (dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, Y] - yq_l_m[m]) * (dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, Y] <= yq_l_m[m])
      Y0.trim.upper <- (dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, Y] - yq_u_m[m]) * (dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, Y] >= yq_u_m[m])
      if (is.null(seed)){
        rf0_l_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, X]), Y = Y0.trim.lower)
        rf0_u_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, X]), Y = Y0.trim.upper)
      }else{
        rf0_l_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, X]), Y = Y0.trim.lower, seed = seed)
        rf0_u_m <- regression_forest(as.matrix(dat.s2[dat.s2[, D] == 0 & dat.s2[, S] == 1, X]), Y = Y0.trim.upper, seed = seed)
      }
      pred0_l_m <- predict(rf0_l_m, newdata = as.matrix(X_moderator), estimate.variance = 1)
      pred0_u_m <- predict(rf0_u_m, newdata = as.matrix(X_moderator), estimate.variance = 1)
      theta0_l_m <- pred0_l_m[m, 1] / one_q_l_m[m] + yq_l_m[m]
      theta0_u_m <- pred0_u_m[m, 1] / one_q_l_m[m] + yq_u_m[m]

      se_theta0_l_m <- sqrt(pred0_l_m[m, 2]) / one_q_l_m[m] + se_one_q_l_m[m] * (pred0_l_m[m, 1] / one_q_l_m[m]^2)
      se_theta0_u_m <- sqrt(pred0_u_m[m, 2]) / one_q_l_m[m] + se_one_q_l_m[m] * (pred0_u_m[m, 1] / one_q_l_m[m]^2)

      tau_l_m[m, k] <- directions_m[m]*(theta1_l_m - theta0_m[m]) + (1-directions_m[m])*(theta1_m[m] - theta0_u_m)
      se_tau_l_m[m, k] <- directions_m[m]*sqrt(se_theta1_l_m^2 + se_theta0_m[m]^2) + (1-directions_m[m])*(se_theta1_m[m]^2 + se_theta0_u_m^2)
      tau_u_m[m, k] <- directions_m[m]*(theta1_u_m - theta0_m[m]) + (1-directions_m[m])*(theta1_m[m] - theta0_l_m)
      se_tau_u_m[m, k] <- directions_m[m]*sqrt(se_theta1_u_m^2 + se_theta0_m[m]^2) + (1-directions_m[m])*(se_theta1_m[m]^2 + se_theta0_l_m^2)
    }
  }


  tau_l_m_avg <- apply(tau_l_m, 1, mean)
  se_tau_l_m_avg <- apply(se_tau_l_m, 1, mean) / sqrt(cv_fold)
  tau_u_m_avg <- apply(tau_u_m, 1, mean)
  se_tau_u_m_avg <- apply(se_tau_l_m, 1, mean) / sqrt(cv_fold)
  tau_m_avg <- apply(tau_m, 1, mean)
  se_tau_m_avg <- apply(se_tau_m, 1, mean) / sqrt(cv_fold)
  result <- list("tau_l_m_avg" = tau_l_m_avg,
                 "tau_u_m_avg" = tau_u_m_avg,
                 "se_tau_l_m_avg" = se_tau_l_m_avg,
                 "se_tau_u_m_avg" = se_tau_u_m_avg)
  return(result)
}

LBounds_cont <- function(dat_final, Y, X, S, D, W, Pscore,
                         cond.mono = FALSE, direction = NULL, message){

  if (cond.mono == 1){
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

    # hurt lee bounds
    all_ys_hurt <- unlist(unique(dat_hurt[dat_hurt[, D] == 0 & dat_hurt[, S] == 1, Y]))
    all_ys_hurt <- sort(all_ys_hurt)

    theta1_lee <- sum(dat_hurt[, D] * dat_hurt[, S] * dat_hurt[, W] * dat_hurt[, Y], na.rm = 1) / sum(dat_hurt[, D] * dat_hurt[, S] * dat_hurt[, W], na.rm = 1)
    var_theta1_lee <- sum((dat_hurt[, D]) * dat_hurt[, S] * dat_hurt[, W] * dat_hurt[, Y]^2, na.rm = 1) / sum(dat_hurt[, D] * dat_hurt[, S] * dat_hurt[, W], na.rm = 1) - theta1_lee^2

    q0_lee_hurt <- sum(dat_hurt[, S] * (1-dat_hurt[, D]) * dat_hurt[, W]) / sum((1-dat_hurt[, D]) * dat_hurt[, W])
    q1_lee_hurt <- sum(dat_hurt[, S] * dat_hurt[, D] * dat_hurt[, W]) / sum(dat_hurt[, D] * dat_hurt[, W])
    q_lee_hurt <- q1_lee_hurt / q0_lee_hurt

    ED_hurt <- pscore_hurt <- sum(dat_hurt[, D] * dat_hurt[, W]) / sum(dat_hurt[, W])

    tilde.F.hurt <- function(y){
      sum(dat_hurt[, S]*(1-dat_hurt[, D])*dat_hurt[, W]*as.numeric(dat_hurt[, Y] < y), na.rm = 1)/sum(dat_hurt[, S]*dat_hurt[, W]*(1-dat_hurt[, D]))
    }

    tilde.F.vals.hurt <- rep(0, length(all_ys_hurt))
    for(i in 1:length(all_ys_hurt)){
      if(i == 1){tilde.F.vals.hurt[i] <- max(0, tilde.F.hurt(all_ys_hurt[i]))}
      if(i > 1){tilde.F.vals.hurt[i] <- max(tilde.F.vals.hurt[i-1], tilde.F.hurt(all_ys_hurt[i]))}
      if(i == length(all_ys_hurt)){tilde.F.vals.hurt[i] <- max(1, tilde.F.hurt(all_ys_hurt[i]))}
    }

    yq_l_lee_hurt <- max(all_ys_hurt[tilde.F.vals.hurt <= q_lee_hurt])
    yq_u_lee_hurt <- min(all_ys_hurt[tilde.F.vals.hurt >= 1-q_lee_hurt])

    theta0_l_lee <- sum((1-dat_hurt[, D])*dat_hurt[, W]*dat_hurt[, S]*as.numeric(dat_hurt[, Y] <= yq_l_lee_hurt)*as.numeric(dat_hurt[, Y]), na.rm = 1) / sum((1-dat_hurt[, D])*dat_hurt[, W]*dat_hurt[, S]*as.numeric(dat_hurt[, Y] <= yq_l_lee_hurt), na.rm = 1)
    theta0_u_lee <- sum((1-dat_hurt[, D])*dat_hurt[, W]*dat_hurt[, S]*as.numeric(dat_hurt[, Y] >= yq_u_lee_hurt)*as.numeric(dat_hurt[, Y]), na.rm = 1) / sum((1-dat_hurt[, D])*dat_hurt[, W]*dat_hurt[, S]*as.numeric(dat_hurt[, Y] >= yq_u_lee_hurt), na.rm = 1)

    tau_u_lee_hurt <- theta1_lee - theta0_l_lee
    tau_l_lee_hurt <- theta1_lee - theta0_u_lee

    V_q_hurt <- q_lee_hurt^2 * (((1 - q1_lee_hurt) / (ED_hurt * q1_lee_hurt)) + (1 - q0_lee_hurt) / ((1-ED_hurt) * q0_lee_hurt))
    ESD_hurt <- sum(dat_hurt[, D] * dat_hurt[, S] * dat_hurt[, W]) / sum(dat_hurt[, W])

    var_tau_u_lee_hurt <- (sum((1-dat_hurt[, D])*dat_hurt[, S]*dat_hurt[, W]*as.numeric(dat_hurt[, Y] <= yq_l_lee_hurt)*as.numeric(dat_hurt[, Y]^2), na.rm = 1) / sum((1-dat_hurt[, D])*dat_hurt[, S]*dat_hurt[, W]*as.numeric(dat_hurt[, Y] <= yq_l_lee_hurt), na.rm = 1) - theta0_l_lee^2) / (ESD_hurt * q_lee_hurt) +
      (1 - q_lee_hurt) * (yq_l_lee_hurt - theta0_l_lee)^2 / (ESD_hurt * q_lee_hurt) + ((yq_l_lee_hurt - theta0_l_lee) / q_lee_hurt)^2 * V_q_hurt + var_theta1_lee
    se_tau_u_lee_hurt <- sqrt(var_tau_u_lee_hurt / nrow(dat_hurt))

    var_tau_l_lee_hurt <- (sum((1-dat_hurt[, D])*dat_hurt[, S]*dat_hurt[, W]*as.numeric(dat_hurt[, Y] >= yq_u_lee_hurt)*as.numeric(dat_hurt[, Y]^2), na.rm = 1) / sum((1-dat_hurt[, D])*dat_hurt[, S]*dat_hurt[, W]*as.numeric(dat_hurt[, Y] >= yq_u_lee_hurt), na.rm = 1) - theta0_u_lee^2) / (ESD_hurt * q_lee_hurt) +
      (1 - q_lee_hurt) * (yq_u_lee_hurt - theta0_u_lee)^2 / (ESD_hurt * q_lee_hurt) + ((yq_u_lee_hurt - theta0_u_lee) / q_lee_hurt)^2 * V_q_hurt + var_theta1_lee
    se_tau_l_lee_hurt <- sqrt(var_tau_l_lee_hurt / nrow(dat_hurt))

    # help lee bounds
    all_ys_help <- unique(dat_help[dat_help[, D] == 1 & dat_help[, S] == 1, Y])
    all_ys_help <- sort(all_ys_help)
    theta0_lee <- sum((1-dat_help[, D]) * dat_help[, S] * dat_help[, W] * dat_help[, Y], na.rm = 1) / sum((1-dat_help[, D]) * dat_help[, S] * dat_help[, W], na.rm = 1)
    var_theta0_lee <- sum((1-dat_help[, D]) * dat_help[, S] * dat_help[, W] * dat_help[, Y]^2, na.rm = 1) / sum((1-dat_help[, D]) * dat_help[, S] * dat_help[, W], na.rm = 1) - theta0_lee^2
    q0_lee_help <- sum(dat_help[, S] * (1-dat_help[, D]) * dat_help[, W]) / sum((1-dat_help[, D]) * dat_help[, W])
    q1_lee_help <- sum(dat_help[, S] * dat_help[, D] * dat_help[, W]) / sum(dat_help[, D] * dat_help[, W])
    q_lee_help <- q0_lee_help / q1_lee_help

    ED_help <- pscore_help <- sum(dat_help[, D] * dat_help[, W]) / sum(dat_help[, W])

    tilde.F.help <- function(y){
      sum(dat_help[, S]*dat_help[, D]*dat_help[, W]*as.numeric(dat_help[, Y] < y), na.rm = 1)/sum(dat_help[, S]*dat_help[, D]*dat_help[, W])
    }

    tilde.F.vals.help <- rep(0, length(all_ys_help))
    for(i in 1:length(all_ys_help)){
      if(i == 1){tilde.F.vals.help[i] <- max(0, tilde.F.help(all_ys_help[i]))}
      if(i > 1){tilde.F.vals.help[i] <- max(tilde.F.vals.help[i-1], tilde.F.help(all_ys_help[i]))}
      if(i == length(all_ys_help)){tilde.F.vals.help[i] <- max(1, tilde.F.help(all_ys_help[i]))}
    }

    yq_l_lee_help <- max(all_ys_help[tilde.F.vals.help <= q_lee_help])
    yq_u_lee_help <- min(all_ys_help[tilde.F.vals.help >= 1-q_lee_help])

    theta1_l_lee <- sum(dat_help[, D]*dat_help[, S]*dat_help[, W]*as.numeric(dat_help[, Y] <= yq_l_lee_help)*as.numeric(dat_help[, Y]), na.rm = 1) / sum(dat_help[, D]*dat_help[, S]*dat_help[, W]*as.numeric(dat_help[, Y] <= yq_l_lee_help), na.rm = 1)
    theta1_u_lee <- sum(dat_help[, D]*dat_help[, S]*dat_help[, W]*as.numeric(dat_help[, Y] >= yq_u_lee_help)*as.numeric(dat_help[, Y]), na.rm = 1) / sum(dat_help[, D]*dat_help[, S]*dat_help[, W]*as.numeric(dat_help[, Y] >= yq_u_lee_help), na.rm = 1)

    tau_u_lee_help <- theta1_u_lee - theta0_lee
    tau_l_lee_help <- theta1_l_lee - theta0_lee

    V_q_help <- q_lee_help^2 * (((1 - q0_lee_help/q_lee_help) / (ED_help * q0_lee_help/q_lee_help)) + (1 - q0_lee_help) / ((1-ED_help) * q0_lee_help))
    ESD_help <- sum(dat_help[, D] * dat_help[, S] * dat_help[, W]) / sum(dat_help[, W])

    var_tau_l_lee_help <- (sum(dat_help[, D]*dat_help[, S]*dat_help[, W]*as.numeric(dat_help[, Y] <= yq_l_lee_help)*as.numeric(dat_help[, Y]^2), na.rm = 1) / sum(dat_help[, D]*dat_help[, S]*dat_help[, W]*as.numeric(dat_help[, Y] <= yq_l_lee_help), na.rm = 1) - theta1_l_lee^2) / (ESD_help * q_lee_help) +
      (1 - q_lee_help) * (yq_l_lee_help - theta1_l_lee)^2 / (ESD_help * q_lee_help) + ((yq_l_lee_help - theta1_l_lee) / q_lee_help)^2 * V_q_help + var_theta0_lee
    se_tau_l_lee_help <- sqrt(var_tau_l_lee_help / nrow(dat_help))

    var_tau_u_lee_help <- (sum(dat_help[, D]*dat_help[, S]*dat_help[, W]*as.numeric(dat_help[, Y] >= yq_u_lee_help)*as.numeric(dat_help[, Y]^2), na.rm = 1) / sum(dat_help[, D]*dat_help[, S]*dat_help[, W]*as.numeric(dat_help[, Y] >= yq_u_lee_help), na.rm = 1) - theta1_u_lee^2) / (ESD_help * q_lee_help) +
      (1 - q_lee_help) * (yq_u_lee_help - theta1_u_lee)^2 / (ESD_help * q_lee_help) + ((yq_u_lee_help - theta1_u_lee) / q_lee_help)^2 * V_q_help + var_theta0_lee
    se_tau_u_lee_help <- sqrt(var_tau_u_lee_help / nrow(dat_help))

    pos_ratio <- mean(directions)

    tau_u_lee <- pos_ratio * tau_u_lee_help + (1-pos_ratio) * tau_u_lee_hurt
    tau_l_lee <- pos_ratio * tau_l_lee_help + (1-pos_ratio) * tau_l_lee_hurt

    se_tau_l_lee <- sqrt(pos_ratio^2 * se_tau_l_lee_help^2 + (1-pos_ratio)^2 * se_tau_l_lee_hurt^2)
    se_tau_u_lee <- sqrt(pos_ratio^2 * se_tau_u_lee_help^2 + (1-pos_ratio)^2 * se_tau_u_lee_hurt^2)
  }else{
    q0_lee <- sum(dat_final[, S] * (1-dat_final[, D]) * dat_final[, W]) / sum((1-dat_final[, D]) * dat_final[, W])
    q1_lee <- sum(dat_final[, S] * dat_final[, D] * dat_final[, W]) / sum(dat_final[, D] * dat_final[, W])

    if (is.null(direction)){
      direction <- as.numeric(q0_lee < q1_lee)
    }
    if (direction == 1){
      if (message == TRUE){
          cat("S(1) > S(0)", "\n")
        }
      all_ys <- dat_final[dat_final[, D] == 1 & dat_final[, S] == 1, Y]
      all_ys_uni <- unique(dat_final[dat_final[, D] == 1 & dat_final[, S] == 1, Y])
      all_ys_uni <- sort(all_ys_uni)
      theta0_lee <- sum((1-dat_final[, D]) * dat_final[, S] * dat_final[, W] * dat_final[, Y], na.rm = 1) / sum((1-dat_final[, D]) * dat_final[, S] * dat_final[, W], na.rm = 1)
      var_theta0_lee <- sum((1-dat_final[, D]) * dat_final[, S] * dat_final[, W] * dat_final[, Y]^2, na.rm = 1) / sum((1-dat_final[, D]) * dat_final[, S] * dat_final[, W], na.rm = 1) - theta0_lee^2

      q_lee <- q0_lee / q1_lee

      ED <- pscore <- sum(dat_final[, D] * dat_final[, W]) / sum(dat_final[, W])

      tilde.F <- function(y){
        sum(dat_final[, S]*dat_final[, D]*dat_final[, W]*as.numeric(dat_final[, Y] < y), na.rm = 1)/sum(dat_final[, S]*dat_final[, D]*dat_final[, W])
      }

      tilde.F.vals <- rep(0, length(all_ys_uni))
      for(i in 1:length(all_ys_uni)){
        if(i == 1){tilde.F.vals[i] <- max(0, tilde.F(all_ys_uni[i]))}
        if(i > 1){tilde.F.vals[i] <- max(tilde.F.vals[i-1], tilde.F(all_ys_uni[i]))}
        if(i == length(all_ys_uni)){tilde.F.vals[i] <- max(1, tilde.F(all_ys_uni[i]))}
      }

      yq_l_lee <- max(all_ys_uni[tilde.F.vals <= q_lee])
      yq_u_lee <- min(all_ys_uni[tilde.F.vals >= 1-q_lee])

      theta1_l_lee <- sum(dat_final[, D]*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] <= yq_l_lee)*as.numeric(dat_final[, Y]), na.rm = 1) / sum(dat_final[, D]*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] <= yq_l_lee), na.rm = 1)
      theta1_u_lee <- sum(dat_final[, D]*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] >= yq_u_lee)*as.numeric(dat_final[, Y]), na.rm = 1) / sum(dat_final[, D]*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] >= yq_u_lee), na.rm = 1)

      tau_u_lee <- theta1_u_lee - theta0_lee
      tau_l_lee <- theta1_l_lee - theta0_lee

      V_q <- q_lee^2 * (((1 - q0_lee/q_lee) / (ED * q0_lee/q_lee)) + (1 - q0_lee) / ((1-ED) * q0_lee))
      ESD <- sum(dat_final[, D] * dat_final[, S] * dat_final[, W]) / sum(dat_final[, W])

      var_tau_l_lee <- (sum(dat_final[, D]*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] <= yq_l_lee)*as.numeric(dat_final[, Y]^2), na.rm = 1) / sum(dat_final[, D]*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] <= yq_l_lee), na.rm = 1) - theta1_l_lee^2) / (ESD * q_lee) +
        (1 - q_lee) * (yq_l_lee - theta1_l_lee)^2 / (ESD * q_lee) + ((yq_l_lee - theta1_l_lee) / q_lee)^2 * V_q + var_theta0_lee
      se_tau_l_lee <- sqrt(var_tau_l_lee / nrow(dat_final))

      var_tau_u_lee <- (sum(dat_final[, D]*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] >= yq_u_lee)*as.numeric(dat_final[, Y]^2), na.rm = 1) / sum(dat_final[, D]*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] >= yq_u_lee), na.rm = 1) - theta1_u_lee^2) / (ESD * q_lee) +
        (1 - q_lee) * (yq_u_lee - theta1_u_lee)^2 / (ESD * q_lee) + ((yq_u_lee - theta1_u_lee) / q_lee)^2 * V_q + var_theta0_lee
      se_tau_u_lee <- sqrt(var_tau_u_lee / nrow(dat_final))
    }else if (direction == 0){
      if (message == TRUE){
        cat("S(0) > S(1)", "\n")
      }
      all_ys <- unlist(dat_final[dat_final[, D] == 0 & dat_final[, S] == 1, Y])
      all_ys_uni <- unlist(unique(dat_final[dat_final[, D] == 0 & dat_final[, S] == 1, Y]))
      all_ys_uni <- sort(all_ys_uni)

      theta1_lee <- sum(dat_final[, D] * dat_final[, S] * dat_final[, W] * dat_final[, Y], na.rm = 1) / sum(dat_final[, D] * dat_final[, S] * dat_final[, W], na.rm = 1)
      var_theta1_lee <- sum((dat_final[, D]) * dat_final[, S] * dat_final[, W] * dat_final[, Y]^2, na.rm = 1) / sum(dat_final[, D] * dat_final[, S] * dat_final[, W], na.rm = 1) - theta1_lee^2

      q_lee <- q1_lee / q0_lee

      ED <- pscore <- sum(dat_final[, D] * dat_final[, W]) / sum(dat_final[, W])

      tilde.F <- function(y){
        sum(dat_final[, S]*(1-dat_final[, D])*dat_final[, W]*as.numeric(dat_final[, Y] < y), na.rm = 1)/sum(dat_final[, S]*(1-dat_final[, D])*dat_final[, W])
      }

      tilde.F.vals <- rep(0, length(all_ys_uni))
      for(i in 1:length(all_ys_uni)){
        if(i == 1){tilde.F.vals[i] <- max(0, tilde.F(all_ys_uni[i]))}
        if(i > 1){tilde.F.vals[i] <- max(tilde.F.vals[i-1], tilde.F(all_ys_uni[i]))}
        if(i == length(all_ys)){tilde.F.vals[i] <- max(1, tilde.F(all_ys_uni[i]))}
      }

      yq_l_lee <- max(all_ys_uni[tilde.F.vals <= q_lee])
      yq_u_lee <- min(all_ys_uni[tilde.F.vals >= 1-q_lee])

      theta0_l_lee <- sum((1-dat_final[, D])*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] <= yq_l_lee)*as.numeric(dat_final[, Y]), na.rm = 1) / sum((1-dat_final[, D])*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] <= yq_l_lee), na.rm = 1)
      theta0_u_lee <- sum((1-dat_final[, D])*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] >= yq_u_lee)*as.numeric(dat_final[, Y]), na.rm = 1) / sum((1-dat_final[, D])*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] >= yq_u_lee), na.rm = 1)

      tau_u_lee <- theta1_lee - theta0_l_lee
      tau_l_lee <- theta1_lee - theta0_u_lee

      V_q <- q_lee^2 * (((1 - q1_lee/q_lee) / ((1-ED) * q1_lee/q_lee)) + (1 - q1_lee) / (ED * q1_lee))
      ESD <- sum((1-dat_final[, D]) * dat_final[, S] * dat_final[, W]) / sum(dat_final[, W])

      var_tau_u_lee <- (sum((1-dat_final[, D])*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] <= yq_l_lee)*as.numeric(dat_final[, Y]^2), na.rm = 1) / sum((1-dat_final[, D])*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] <= yq_l_lee), na.rm = 1) - theta0_l_lee^2) / (ESD * q_lee) +
        (1 - q_lee) * (yq_l_lee - theta0_l_lee)^2 / (ESD * q_lee) + ((yq_l_lee - theta0_l_lee) / q_lee)^2 * V_q + var_theta1_lee
      se_tau_u_lee <- sqrt(var_tau_u_lee / nrow(dat_final))

      var_tau_l_lee <- (sum((1-dat_final[, D])*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] >= yq_u_lee)*as.numeric(dat_final[, Y]^2), na.rm = 1) / sum((1-dat_final[, D])*dat_final[, S]*dat_final[, W]*as.numeric(dat_final[, Y] >= yq_u_lee), na.rm = 1) - theta0_u_lee^2) / (ESD * q_lee) +
        (1 - q_lee) * (yq_u_lee - theta0_u_lee)^2 / (ESD * q_lee) + ((yq_u_lee - theta0_u_lee) / q_lee)^2 * V_q + var_theta1_lee
      se_tau_l_lee <- sqrt(var_tau_l_lee / nrow(dat_final))
    }else{
      stop("Direction can only be 0 or 1")
    }
  }
  result <- list("tau_l_est_lee" = tau_l_lee,
                 "tau_u_est_lee" = tau_u_lee,
                 "se_tau_l_est_lee" = se_tau_l_lee,
                 "se_tau_u_est_lee" = se_tau_u_lee)
  return(result)
}

