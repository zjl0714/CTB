CTTB <- function(data, seed = NULL, Y, D, S, X = NULL, W = NULL, Pscore = NULL,
                 regression.splitting = FALSE, cv_fold = 5, trim_l = 0, trim_u = 1,
                 aggBounds = TRUE, IM_cv = TRUE, alpha = 0.05, cBounds = FALSE, X_moderator = NULL, 
                 LeeBounds = FALSE, direction = NULL, cond.mono = TRUE, message = TRUE){
  
  if (is.null(data)){
    stop("No dataset specified")
  }
  
  if (!is.null(seed)){
    set.seed(seed)
  }

  if (is.null(D)){
    stop("No treatment specified")
  }else if (!all(data[, D] %in% 0:1)){
    stop("Non-binary treatment")
  }
  
  if (is.null(Pscore) & message == TRUE){
    cat("Propensity scores are not specified and will be estimated using X.", "\n")
  }
  
  if (is.null(W)){
    W <- "uni_w"
    data[ W] <- 1
  }
  N <- sum(data[, W])
  
  if (is.null(S)){
    stop("No missing/selection indicator specified")
  }else if (!all(data[, S] %in% 0:1)){
    stop("The missing/selection indicator is not binary")
  }else{
    S1_sum <- sum(data[, S])
    S0_sum <- sum(1-data[, S])
    if (message == TRUE){
      cat("The average missing rate is ", round(S0_sum/N, 4), "\n")
    }
  }
  
  bo <- 0
  if (is.null(Y)){
    stop("No outcome specified")
  }else if (all(data[data[, S] == 1, Y] %in% 0:1)){
    cat("Binary outcome", "\n")
    bo <- 1
  }
  
  
  if (is.null(X)){
    cat("No covariate. Lee bounds will be used", "\n")
    aggBounds <- FALSE
    LeeBounds <- TRUE
  }else{
    P <- dim(data[, X])[2]
    for (pp in 1:P){
      X_p <- data[, X[pp]]
      if (!is.numeric(X_p)){
        stop("Non-numeric covariates")
      }
    }
  }
  
  data0 <- data[data[, D] == 0, ]
  data1 <- data[data[, D] == 1, ]
  
  folds_0 <- cut(seq(1, nrow(data0)), breaks = cv_fold, labels = FALSE)
  folds_1 <- cut(seq(1, nrow(data1)), breaks = cv_fold, labels = FALSE)
  
  
  # check splitting and ensure there are more than 
  # 5 treated and untreated obs in each fold
  flag <- TRUE
  while (flag){
    split_index_0 <- sample(folds_0, length(folds_0), replace = F)
    split_index_1 <- sample(folds_1, length(folds_1), replace = F)
    
    data0$split_index <- split_index_0
    data1$split_index <- split_index_1
    
    dat_final <- data.frame(rbind(data0, data1))
    flags <- rep(NA, cv_fold)
    for (cc in 1:cv_fold){
      d_train <- dat_final[dat_final$split_index != cc, ]
      flags[cc] <- sum(d_train[d_train[, D] == 1, S]) < 5 | 
        sum(1 - d_train[d_train[, D] == 1, S]) < 5 | 
        sum(d_train[d_train[, D] == 0, S]) < 5 | 
        sum(1 - d_train[d_train[, D] == 0, S]) < 5
    }
    flag <- mean(flags) > 0
  }
  
  result <- list()
  if (bo == 1){
    if (LeeBounds){
      lee.result <- LBounds_binary(dat_final, Y, X, S, D, W, Pscore,
                                   cond.mono, direction)
      if (IM_cv){
        IM_cv_lee <- IM_cv(lee.result[["tau_l_est_lee"]], lee.result[["tau_u_est_lee"]], 
                           lee.result[["se_tau_l_est_lee"]], lee.result[["se_tau_u_est_lee"]], 
                           alpha, sum(dat[, S]))
      }
      TB_dat <- data.frame("Estimate" = c(lee.result[["tau_l_est_lee"]], lee.result[["tau_u_est_lee"]]),
                           "SE" = c(lee.result[["se_tau_l_est_lee"]], lee.result[["se_tau_u_est_lee"]]))
      TB_dat$CI95_lower <- TB_dat$Estimate - qnorm(1 - alpha/2) * TB_dat$SE
      TB_dat$CI95_upper <- TB_dat$Estimate + qnorm(1 - alpha/2) * TB_dat$SE
      if (IM_cv){
        TB_dat$CI95_lower_IM <- TB_dat$Estimate - IM_cv_lee * TB_dat$SE
        TB_dat$CI95_upper_IM <- TB_dat$Estimate + IM_cv_lee * TB_dat$SE
      }
      TB_dat$Method <- rep(c("TB"), 2)
      TB_dat$Type <- c("Lower", "Upper")
      
      result[["TB"]] <- TB_dat
    }
    
    if (aggBounds){
      agg.result <- aggBounds_binary(dat_final = dat_final, seed = seed, Y, X, S, D, W, Pscore,
                                     regression.splitting = regression.splitting,
                                     cv_fold = cv_fold, trim_l, trim_u)
      if (IM_cv){
        IM_cv_agg <- IM_cv(agg.result[["tau_l_est_avg"]], agg.result[["tau_u_est_avg"]], 
                           agg.result[["se_tau_l_est_avg"]], agg.result[["se_tau_u_est_avg"]], 
                           alpha, sum(dat[, S]))
      }
      CTTB_dat <- data.frame("Estimate" = c(agg.result[["tau_l_est_avg"]], agg.result[["tau_u_est_avg"]]),
                             "SE" = c(agg.result[["se_tau_l_est_avg"]], agg.result[["se_tau_u_est_avg"]]))
      CTTB_dat$CI95_lower <- CTTB_dat$Estimate - qnorm(1 - alpha/2) * CTTB_dat$SE
      CTTB_dat$CI95_upper <- CTTB_dat$Estimate + qnorm(1 - alpha/2) * CTTB_dat$SE
      if (IM_cv){
        CTTB_dat$CI95_lower_IM <- CTTB_dat$Estimate - IM_cv_lee * CTTB_dat$SE
        CTTB_dat$CI95_upper_IM <- CTTB_dat$Estimate + IM_cv_lee * CTTB_dat$SE
      }
      CTTB_dat$Method <- rep(c("CTTB"), 2)
      CTTB_dat$Type <- c("Lower", "Upper")
      
      result[["CTTB"]] <- CTTB_dat
    }
  }else{
    if (LeeBounds){
      lee.result <- LBounds_cont(dat_final, Y, X, S, D, W, Pscore,
                                 cond.mono, direction, message=message)
      if (IM_cv){
        IM_cv_lee <- IM_cv(lee.result[["tau_l_est_lee"]], lee.result[["tau_u_est_lee"]], 
                           lee.result[["se_tau_l_est_lee"]], lee.result[["se_tau_u_est_lee"]], 
                           alpha, sum(dat[, S]))
      }
      TB_dat <- data.frame("Estimate" = c(lee.result[["tau_l_est_lee"]], lee.result[["tau_u_est_lee"]]),
                           "SE" = c(lee.result[["se_tau_l_est_lee"]], lee.result[["se_tau_u_est_lee"]]))
      TB_dat$CI95_lower <- TB_dat$Estimate - qnorm(1 - alpha/2) * TB_dat$SE
      TB_dat$CI95_upper <- TB_dat$Estimate + qnorm(1 - alpha/2) * TB_dat$SE
      if (IM_cv){
        TB_dat$CI95_lower_IM <- TB_dat$Estimate - IM_cv_lee * TB_dat$SE
        TB_dat$CI95_upper_IM <- TB_dat$Estimate + IM_cv_lee * TB_dat$SE 
      }
      TB_dat$Method <- rep(c("TB"), 2)
      TB_dat$Type <- c("Lower", "Upper")
      
      result[["TB"]] <- TB_dat
    }
    
    if (aggBounds){
      agg.result <- aggBounds_cont(dat_final, seed = seed, Y, X, S, D, W, Pscore,
                                   regression.splitting = regression.splitting,
                                   cv_fold = cv_fold, trim_l, trim_u)
      if (IM_cv){
        IM_cv_agg <- IM_cv(agg.result[["tau_l_est_avg"]], agg.result[["tau_u_est_avg"]], 
                           agg.result[["se_tau_l_est_avg"]], agg.result[["se_tau_u_est_avg"]], 
                           alpha, sum(dat[, S]))
      }
      CTTB_dat <- data.frame("Estimate" = c(agg.result[["tau_l_est_avg"]], agg.result[["tau_u_est_avg"]]),
                             "SE" = c(agg.result[["se_tau_l_est_avg"]], agg.result[["se_tau_u_est_avg"]]))
      CTTB_dat$CI95_lower <- CTTB_dat$Estimate - qnorm(1 - alpha/2) * CTTB_dat$SE
      CTTB_dat$CI95_upper <- CTTB_dat$Estimate + qnorm(1 - alpha/2) * CTTB_dat$SE
      if (IM_cv){
        CTTB_dat$CI95_lower_IM <- CTTB_dat$Estimate - IM_cv_lee * CTTB_dat$SE
        CTTB_dat$CI95_upper_IM <- CTTB_dat$Estimate + IM_cv_lee * CTTB_dat$SE 
      }
      CTTB_dat$Method <- rep(c("CTTB"), 2)
      CTTB_dat$Type <- c("Lower", "Upper")
      
      result[["CTTB"]] <- CTTB_dat
      result[["Ps_k"]] <- agg.result[["Ps_k"]]
    }
    
    if (cBounds){
      if (is.null(X_moderator)){
        stop("No moderator is provided.")
      }
      cond.result <- condBounds_cont(dat_final, seed = seed, Y, X, S, D, W, Pscore,
                                     X_moderator, regression.splitting = regression.splitting,
                                     cv_fold = cv_fold)
      # if (IM_cv){
      #   IM_cv_cond <- rep(NA, nrow(X_moderator))
      #   for (rr in 1:nrow(X_moderator)){
      #     IM_cv_cond[rr] <- IM_cv(cond.result[["tau_l_est_avg"]][rr], cond.result[["tau_u_est_avg"]][rr], 
      #                             cond.result[["se_tau_l_est_avg"]][rr], cond.result[["se_tau_u_est_avg"]][rr], 
      #                             alpha, sum(dat[, S]))
      #   }
      # }
      cond_dat <- data.frame("Estimate" = c(cond.result[["tau_l_m_avg"]], cond.result[["tau_u_m_avg"]]),
                             "SE" = c(cond.result[["se_tau_l_m_avg"]], cond.result[["se_tau_u_m_avg"]]))
      cond_dat$CI95_lower <- cond_dat$Estimate - qnorm(1 - alpha/2) * cond_dat$SE
      cond_dat$CI95_upper <- cond_dat$Estimate + qnorm(1 - alpha/2) * cond_dat$SE
      # if (IM_cv){
      #   cond_dat$CI95_lower_IM <- cond_dat$Estimate - IM_cv_cond * cond_dat$SE
      #   cond_dat$CI95_upper_IM <- cond_dat$Estimate + IM_cv_cond * cond_dat$SE
      # }
      cond_dat$Method <- rep(c("Conditional Bounds"), 2*nrow(X_moderator))
      cond_dat$Type <- rep(c("Lower", "Upper"), each = nrow(X_moderator))
      colnames(X_moderator) <- c("Const", paste0("X", seq(1:(ncol(X_moderator) - 1))))
      cond_dat <- cbind(cond_dat, X_moderator)
      
      result[["cBounds"]] <- cond_dat
    }
  }
  result[["parameters"]] <- c("bo" = bo, "LeeBounds" = LeeBounds, "aggBounds" = aggBounds, 
                              "cBounds" = cBounds, "IM_cv" = IM_cv, "cond.mono" = cond.mono)
  class(result) <- "CTTB"
  result$call <- match.call()
  return(result)
}
 
 
IM_cv <- function(tau_l, tau_u, se_tau_l, se_tau_u, alpha, N){
  cv.fun <- function(cv, tau_l, tau_u, se_tau_l, se_tau_u, alpha, N){
    pnorm(cv + (tau_u - tau_l) / max(se_tau_l, se_tau_u)) - pnorm(-cv) - (1 - alpha)
  }
  cv.IM <- try(uniroot(cv.fun, interval=c(-5, 5), tau_l, tau_u, se_tau_l, se_tau_u, alpha, N), silent = TRUE)
  return(cv.IM$root)
}

print.CTTB <- function(result, ...) {
  if (result[["parameters"]]["LeeBounds"]){
    cat("Estimates from trimming bounds:", "\n")
    print(result[["TB"]])
  }
  if (result[["parameters"]]["aggBounds"]){
    cat("Estimates from covariate-tightened trimming bounds:", "\n")
    print(result[["CTTB"]])
  }
  if (result[["parameters"]]["cBounds"]){
    cat("Estimates of conditional bounds:", "\n")
    print(result[["cBounds"]])
  }
  invisible(result)
}

summary.CTTB <- function(result, ...) {
  if (result[["parameters"]]["LeeBounds"]){
    cat("Estimates from trimming bounds:", "\n")
    print(result[["TB"]])
  }
  if (result[["parameters"]]["aggBounds"]){
    cat("Estimates from covariate-tightened trimming bounds:", "\n")
    print(result[["CTTB"]])
  }
  if (result[["parameters"]]["cBounds"]){
    cat("Estimates of conditional bounds:", "\n")
    print(result[["cBounds"]])
  }
  invisible(result)
}

plot.CTTB <- function(result, type = "aggregated", moderator_plot = 1, ...) {
  if (type == "aggregated"){
    if (result[["parameters"]]["LeeBounds"]){
      p.data <- result[["TB"]]
      if (result[["parameters"]]["aggBounds"]){
        p.data <- rbind(p.data, result[["CTTB"]])
      }
    }else if (result[["parameters"]]["aggBounds"]){
      p.data <- result[["CTTB"]]
    }
    pd <- position_dodge(width=0.2)
    coefs_p <- ggplot(p.data, aes(Method, Estimate, color=Type)) +
      geom_point(aes(shape=Method), size=4, position=pd, show.legend = FALSE) +
      scale_color_manual(name="Estimand",values=c("grey", "black")) +
      scale_shape_manual(name="Method",values=c(19, 19)) +
      theme_bw() +
      geom_errorbar(aes(Method, ymin=CI95_lower, ymax=CI95_upper), width=0.2, position=pd, linetype = "dashed") + 
      xlab('') + ylab('') +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), axis.text.x = element_text(angle = 0), legend.position = "top") + 
      ylab("Estimates") + coord_flip()
    
    y_min <- min(p.data$CI95_upper, na.rm = TRUE)
    
    if (result[["parameters"]]["IM_cv"]){
      coefs_p <- coefs_p + geom_errorbar(aes(Method, ymin=CI95_lower_IM, ymax=CI95_upper_IM), width=0.2, position=pd) +
        annotate("segment", x = 2.2, xend = 2.2, y = y_min - 0.6, yend = y_min - 0.2, linetype = "dashed", size = 0.5) +
        annotate("text", x = 2.2, y = y_min, label = "CI for the Bounds", hjust = 0, size = 5) +
        annotate("segment", x = 2.5, xend = 2.5, y = y_min - 0.6, yend = y_min - 0.2, linetype = "solid", size = 0.5) +
        annotate("text", x = 2.5, y = y_min, label = "IM CI for the Estimate", hjust = 0, size = 5)
    }else{
      coefs_p <- coefs_p + 
        annotate("segment", x = 2.2, xend = 2.2, y = y_min - 0.6, yend = y_min - 0.2, linetype = "dashed", size = 0.5) +
        annotate("text", x = 2.2, y = y_min, label = "CI for the Bounds", hjust = 0, size = 5)
    }
    return(coefs_p)
  }else if (type == "conditional"){
    if (result[["parameters"]]["cBounds"]){
      p.data.m <- result[["cBounds"]]
      p.data.m$Type <- factor(p.data.m$Type, levels = c("Upper", "Lower"),
                              labels = c("Upper Bound", "Lower Bound"))
      moderator_name <- sym(paste0("X", moderator_plot))
      if (var(p.data.m[, paste0("X", moderator_plot)]) == 0){
        stop("No variation in the moderator. Plots cannot be produced.")
      }
      coefs_p <- ggplot(p.data.m, aes(x = !!moderator_name, y = Estimate, group = Type)) +
        geom_line(aes(y = CI95_lower, linetype = Type, color = Type), size = 1) +
        geom_line(aes(y = CI95_upper, linetype = Type, color = Type), size = 1) +
        geom_point(aes(color = Type), size = 3) +
        labs(x = "Moderator", y = "Conditional Bounds", color = NULL, linetype = NULL) +
        scale_color_manual(values = c("black", "gray")) +
        scale_linetype_manual(values = c("dashed", "dashed")) +
        theme_bw(base_size = 16) +
        theme(legend.position = "top")
    }
    return(coefs_p)
  }
}
