CTTB <- function(data, seed = NULL, Y, D, S, X = NULL, W = NULL, Pscore = NULL,
                 regression.splitting = FALSE, cv_fold = 5, trim_l = 0, trim_u = 1,
                 aggBounds = TRUE,
                 cBounds = FALSE, X_moderator = NULL, 
                 LeeBounds = FALSE, direction = NULL,
                 cond.mono = TRUE){
  
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
  
  if (is.null(Pscore)){
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
    cat("The average missing rate is ", round(S0_sum/N, 4), "\n")
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
    result <- LeeBounds_cont()
    return(result)
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
  
  folds_0 <- cut(seq(1, nrow(data0)), breaks=cv_fold, labels=FALSE)
  folds_1 <- cut(seq(1, nrow(data1)), breaks=cv_fold, labels=FALSE)
  
  split_index_0 <- sample(folds_0, length(folds_0), replace = F)
  split_index_1 <- sample(folds_1, length(folds_1), replace = F)
  
  data0$split_index <- split_index_0
  data1$split_index <- split_index_1
  
  dat_final <- data.frame(rbind(data0, data1))
  
  result <- list()
  if (bo == 1){
    if (LeeBounds){
      lee.result <- LBounds_binary(dat_final, Y, X, S, D, W, Pscore,
                                   cond.mono, direction)
      result[["tau_l_est_lee"]] <- lee.result[["tau_l_est_lee"]]
      result[["tau_u_est_lee"]] <- lee.result[["tau_u_est_lee"]]
      result[["se_tau_l_est_lee"]] <- lee.result[["se_tau_l_est_lee"]]
      result[["se_tau_u_est_lee"]] <- lee.result[["se_tau_u_est_lee"]]
    }
    
    if (aggBounds){
      agg.result <- aggBounds_binary(dat_final = dat_final, seed = seed, Y, X, S, D, W, Pscore,
                                     regression.splitting = regression.splitting,
                                     cv_fold = cv_fold, trim_l, trim_u)
      
      result[["tau_l_est_avg"]] <- agg.result[["tau_l_est_avg"]]
      result[["tau_u_est_avg"]] <- agg.result[["tau_u_est_avg"]]
      
      result[["se_tau_l_est_avg"]] <- agg.result[["se_tau_l_est_avg"]]
      result[["se_tau_u_est_avg"]] <- agg.result[["se_tau_u_est_avg"]]
    }
  }else{
    if (LeeBounds){
      lee.result <- LBounds_cont(dat_final, Y, X, S, D, W, Pscore,
                                 cond.mono, direction)
      result[["tau_l_est_lee"]] <- lee.result[["tau_l_est_lee"]]
      result[["tau_u_est_lee"]] <- lee.result[["tau_u_est_lee"]]
      result[["se_tau_l_est_lee"]] <- lee.result[["se_tau_l_est_lee"]]
      result[["se_tau_u_est_lee"]] <- lee.result[["se_tau_u_est_lee"]]
    }
    
    if (aggBounds){
      agg.result <- aggBounds_cont(dat_final, seed = seed, Y, X, S, D, W, Pscore,
                                   regression.splitting = regression.splitting,
                                   cv_fold = cv_fold, trim_l, trim_u)
      
      result[["tau_l_est_avg"]] <- agg.result[["tau_l_est_avg"]]
      result[["tau_u_est_avg"]] <- agg.result[["tau_u_est_avg"]]
      
      result[["se_tau_l_est_avg"]] <- agg.result[["se_tau_l_est_avg"]]
      result[["se_tau_u_est_avg"]] <- agg.result[["se_tau_u_est_avg"]]
      
      result[["Ps_k"]] <- agg.result[["Ps_k"]]
    }
    
    if (cBounds){
      if (is.null(X_moderator)){
        X_moderator <- as.matrix(dat_final[, X])
      }
      cond.result <- condBounds_cont(dat_final, seed = seed, Y, X, S, D, W, Pscore,
                                     X_moderator, regression.splitting = regression.splitting,
                                     cv_fold = cv_fold)
      result[["tau_l_m_avg"]] <- cond.result[["tau_l_m_avg"]]
      result[["tau_u_m_avg"]] <- cond.result[["tau_u_m_avg"]]
      
      result[["se_tau_l_m_avg"]] <- cond.result[["se_tau_l_m_avg"]]
      result[["se_tau_u_m_avg"]] <- cond.result[["se_tau_u_m_avg"]]
    }
  }
  return(result)
}
 
 


