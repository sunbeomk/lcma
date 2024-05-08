# variable selection
library(mclust)




var_select <- function(data, B = 1000, cr = "TIE") {
  
  
  Y <- data$Y
  M <- data$M
  G <- data$G
  K <- ncol(M)
  
  
  # set the initial variable set ----
  # initial clustering with the full data
  mclust_fit_full <- mclust::Mclust(data = M, G = 2:6, verbose = FALSE, modelNames = "VEI")
  full_L <- mclust_fit_full$G
  prob_var <- rep(NA, K)
  for (k in 1:K) {
    mclust_fit_single <- mclust::Mclust(data = M[, k], G = full_L, verbose = FALSE, modelNames = "V")
    if (is.null(mclust_fit_single)) {
      prob_var[k] <- 0
    } else { 
      prob_var[k] <- mean(apply(mclust_fit_single$z, 1, sd))
      }
  }
  selected_list <- remaining_list <- list()
  selected_list[[1]] <- selected <- order(prob_var, decreasing = TRUE)[1:(full_L - 1)]
  remaining_list[[1]] <- remaining <- (1:K)[-selected]
  boot_results <- list()
  
  
  
  
  # clustering with initial variable set
  if(length(selected) > 1) {model_initial <- "VEI"} else {model_initial <- NULL}
  mclust_fit_initial <- mclust::Mclust(data = M[, selected], G = 2:6, verbose = FALSE, modelNames = model_initial)
  Z0_no <- mclust_fit_initial$classification
  Z0_temp <- matrix(0, nrow = N, ncol = mclust_fit_initial$G)
  for (i in 1:N) {
    Z0_temp[i, Z0_no[i]] <- 1
  }
  Z0_initial <- Z0_temp
  if (0 %in% colSums(Z0_initial)) {
    Z0_initial <- Z0_initial[, -which(colSums(Z0_initial) == 0)]
  }
  L_initial <- ncol(Z0_initial)
  # fit lcma model with initial variable set
  fit_initial <- lcma(Y = Y, M = M[, selected], G = G, L = L_initial, Z0_initial)
  
  TIEs <- fit_initial$est$TIE
  DEs <- fit_initial$est$DE
  # bootstrap CI
  if (cr %in% c("TIE")){
    boot_result <- boot_results[[1]] <- bootstrap(fit_initial, B = B)
  }
  if (cr %in% c("em_TIE")){
    boot_result <- boot_results[[1]] <- emp_bootstrap(fit_initial, B = B)
  }
  TIE_CIs <- TIE_CI <- t(quantile(boot_result[, "TIE"], c(0.025, 0.975), na.rm = TRUE))
  DE_CIs <- DE_CI <- t(quantile(boot_result[, "DE"], c(0.025, 0.975), na.rm = TRUE))
  
  
  
  
  # start variable selection ----
  iter <- 2
  while (TRUE) {
    inclusion <- exclusion <- FALSE
    TIE_CIs <- rbind(TIE_CIs, c(NA, NA))
    DE_CIs <- rbind(DE_CIs, c(NA, NA))
    
    
    ## inclusion step ----
    TIE_vec <- DE_vec <- rep(NA, length(remaining))
    for (r in 1:length(remaining)) {
      
      # subset of M including the candidate variable
      select_temp <- c(selected, remaining[r])
      M_temp <- M[, select_temp]
      
      # initial values from mclust
      if(length(select_temp) > 1) {model_temp <- "VEI"} else {model_temp <- NULL}
      mclust_fit_temp <- mclust::Mclust(data = M_temp, G = 2:6, verbose = FALSE, modelNames = model_temp)
      Z0_no <- mclust_fit_temp$classification
      Z0_temp <- matrix(0, nrow = N, ncol = mclust_fit_temp$G)
      for (i in 1:N) {
        Z0_temp[i, Z0_no[i]] <- 1
      }
      Z0 <- Z0_temp
      
      # fit lcma
      tryCatch({fit_temp <- lcma(Y = Y, M = M_temp, G = G, L = mclust_fit_temp$G, Z0)
      TIE_vec[r] <- fit_temp$est$TIE
      DE_vec[r] <- fit_temp$est$DE},
      error = function(e) 0)
      
    }
    TIEs[iter] <- TIE_vec[which.max(abs(TIE_vec))]
    DEs[iter] <- DE_vec[which.min(abs(DE_vec))]
    
    
    # check if absolute value of TIE increases after adding an indicator
    if (abs(TIEs[iter - 1]) < max(abs(TIE_vec), na.rm = TRUE)) { 
      add_ind <- remaining[which.max(abs(TIE_vec))]
      select_temp <- c(selected, add_ind)
      
      # fit the adjusted model
      M_max <- M[, select_temp]
      if(length(select_temp) > 1) {model_max <- "VEI"} else {model_max <- NULL}
      mclust_fit_max <- mclust::Mclust(data = M_max, G = 2:6, verbose = FALSE, modelNames = model_max)
      Z0_no <- mclust_fit_max$classification
      Z0_temp <- matrix(0, nrow = N, ncol = mclust_fit_max$G)
      for (i in 1:N) {
        Z0_temp[i, Z0_no[i]] <- 1
      }
      Z0 <- Z0_temp
      fit_max <- lcma(Y = Y, M = M_max, G = G, L = mclust_fit_max$G, Z0)
      
      # Bootstrap quantile confidence interval
      if (cr %in% c("TIE")) {
        boot_result <- boot_results[[iter]] <- bootstrap(fit_max, B = B)
      }
      if (cr %in% c("em_TIE")) {
        boot_result <- boot_results[[iter]] <- emp_bootstrap(fit_max, B = B)
      }
      TIE_CI <- quantile(boot_result[, "TIE"], c(0.025, 0.975), na.rm = TRUE)
      DE_CI <- quantile(boot_result[, "DE"], c(0.025, 0.975), na.rm = TRUE)
      
      # TIE / DE difference test
      TIE_diff_p <- 1 - mean(abs(boot_result[, "TIE"]) > abs(boot_results[[iter - 1]][, "TIE"]), na.rm = TRUE)
      DE_diff_p <- 1 - mean(abs(boot_result[, "DE"]) < abs(boot_results[[iter - 1]][, "DE"]), na.rm = TRUE)
      
      if (TIE_diff_p < 0.05) { # if the TIE difference is significant
        inclusion <- TRUE
        
        selected_list[[iter]] <- selected <- c(selected, add_ind)
        remaining_list[[iter]] <- remaining <- remaining[-which.max(abs(TIE_vec))]
        
        TIE_CIs[iter, ] <- TIE_CI
        DE_CIs[iter, ] <- DE_CI
      } 
      
    }
    
    if (!inclusion) {
      TIEs[iter] <- TIEs[iter - 1]
      boot_results[[iter]] <- boot_results[[iter - 1]]
    }
    
    
    if (!inclusion & length(selected) == 1) {break}  # break if a single variable is optimal
    
    
    
    
    ## exclusion step ----
    TIE_vec <- DE_vec <- rep(NA, length(selected))
    for (r in 1:length(selected)) {
      
      # subset of M excluding the candidate variable
      select_temp <- selected[-r]
      M_temp <- M[, select_temp]
      
      # initial values from mclust
      if(length(select_temp) > 1) {model_temp <- "VEI"} else {model_temp <- NULL}
      mclust_fit_temp <- mclust::Mclust(data = M_temp, G = 2:6, verbose = FALSE, modelNames = model_temp)
      Z0_no <- mclust_fit_temp$classification
      Z0_temp <- matrix(0, nrow = N, ncol = mclust_fit_temp$G)
      for (i in 1:N) {
        Z0_temp[i, Z0_no[i]] <- 1
      }
      Z0 <- Z0_temp
        
      # fit lcma
      tryCatch({fit_temp <- lcma(Y = Y, M = M_temp, G = G, L = mclust_fit_temp$G, Z0)
      TIE_vec[r] <- fit_temp$est$TIE
      DE_vec[r] <- fit_temp$est$DE},
      error = function(e) 0)

    }
    
    
    # check if absolute value of TIE increases after adding in indicator
    if (abs(TIEs[iter]) < max(abs(TIE_vec), na.rm = TRUE)) {
      remove_ind <- selected[which.max(abs(TIE_vec))]
      select_temp <- selected[-which.max(abs(TIE_vec))]
      
      # fit the adjusted model
      M_max <- M[, select_temp]
      if(length(select_temp) > 1) {model_max <- "VEI"} else {model_max <- NULL}
      mclust_fit_max <- mclust::Mclust(data = M_max, G = 2:6, verbose = FALSE, modelNames = model_max)
      Z0_no <- mclust_fit_max$classification
      Z0_temp <- matrix(0, nrow = N, ncol = mclust_fit_max$G)
      for (i in 1:N) {
        Z0_temp[i, Z0_no[i]] <- 1
      }
      Z0 <- Z0_temp
      if (0 %in% colSums(Z0)) {Z0 <- Z0[, -which(colSums(Z0) == 0)]}
      fit_max <- lcma(Y = Y, M = M_max, G = G, L = ncol(Z0), Z0)
      
      # Bootstrap quantile confidence interval
      if (cr %in% c("TIE")) {
        boot_result <- bootstrap(fit_max, B = B)
      }
      if (cr %in% c("em_TIE")) {
        boot_result <- emp_bootstrap(fit_max, B = B)
      }
      TIE_CI <- quantile(boot_result[, "TIE"], c(0.025, 0.975), na.rm = TRUE)
      DE_CI <- quantile(boot_result[, "DE"], c(0.025, 0.975), na.rm = TRUE)
      
      # TIE / DE difference test
      TIE_diff_p <- 1 - mean(abs(boot_result[, "TIE"]) > abs(boot_results[[iter]][, "TIE"]), na.rm = TRUE)
      DE_diff_p <- 1 - mean(abs(boot_result[, "DE"]) < abs(boot_results[[iter]][, "DE"]), na.rm = TRUE)
      
      if (TIE_diff_p < 0.05) { # if the TIE difference is significant
        exclusion <- TRUE
        
        selected_list[[iter]] <- selected <- selected[-which.max(abs(TIE_vec))]
        remaining_list[[iter]] <- remaining <- c(remaining, remove_ind)
        
        TIE_CIs[iter, ] <- TIE_CI
        DE_CIs[iter, ] <- DE_CI
      } 
      
    }

    
    
    
    print(iter)
    if (!inclusion & !exclusion) {break}
    iter <- iter + 1
    if (length(selected) == K) {break}
    
    
  } # end of variable selection algorithm
  
  
  
  
  # Fit the final model ----
  M_final <- data$M[, selected_list[[iter - 1]]]
  if(length(selected_list[[iter - 1]]) > 1) {model_final <- "VEI"} else {model_final <- NULL}
  mclust_fit_final <- mclust::Mclust(data = M_final, G = 2:6, verbose = FALSE, modelNames = model_final)
  Z0_no <- mclust_fit_final$classification
  Z0_temp <- matrix(0, nrow = N, ncol = mclust_fit_final$G)
  for (i in 1:N) {
    Z0_temp[i, Z0_no[i]] <- 1
  }
  Z0_final <- Z0_temp
  if (0 %in% colSums(Z0_final)) {
    Z0_final <- Z0_final[, -which(colSums(Z0_final) == 0)]
  }
  L_final <- ncol(Z0_final)
  fit_final <- lcma(Y = Y, M = M_final, G = G, L = L_final, Z0_final)
  
  
  
  
  # return results ----
  list(selected = selected_list[1:(iter - 1)],
       iter = iter - 1,
       TIEs = TIEs[1:(iter - 1)],
       DEs = DEs[1:(iter - 1)],
       fit_final = fit_final,
       cr = cr)
  
  
}




