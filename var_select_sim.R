## Load packages ----
library(combinat)




## source functions ----
src_path <- file.path("/home", "sunbeom2", "scratch", "lcma")
source(file.path(src_path, "functions.R"))
source(file.path(src_path, "LCMA_data_gen.R"))
source(file.path(src_path, "LCMA_EM.R"))
source(file.path(src_path, "var_select.R"))
args = commandArgs(trailingOnly = TRUE)

# source("functions.R")
# source("LCMA_data_gen.R")
# source("LCMA_EM.R")
# source("var_select.R")



## Simulation Setting ----------------------------------------------------------
reps <- 2            # number of replications
Ns <- c(500)         # sample size
Ls <- c(4)           # number of latent classes
Ks <- c(10)          # number of process features (M)
NNs <- c(5, 7, 9)    # number of noise variables
VARs <- c(1, 3)      # class variances size proportions


condition <- expand.grid(N = Ns, L = Ls, K = Ks, NN = NNs, VAR = VARs)


L <- 4
K <- 10
true_mus_list <- list()
NN <- 5
true_mus_list[[1]] <- true_mus_list[[4]] <- rbind(
  cbind(rep(0, L - 1), 10 * diag(L - 1)),
  matrix(c(10, rep(0, L - 1)), nrow = K - NN - L + 1, ncol = L, byrow = TRUE),
  matrix(-10, nrow = NN, ncol = L)
  )
NN <- 7
true_mus_list[[2]] <- true_mus_list[[5]] <- rbind(
  cbind(rep(0, L - 1), 10 * diag(L - 1)),
  matrix(c(10, rep(0, L - 1)), nrow = K - NN - L + 1, ncol = L, byrow = TRUE),
  matrix(-10, nrow = NN, ncol = L)
)
NN <- 9
true_mus_list[[3]] <- true_mus_list[[6]] <- rbind(10 * (1:4), matrix(-10, nrow = NN, ncol = L))




## run simulation --------------------------------------------------------------
for (cc in 1:nrow(condition)) {
  
  start <- Sys.time()
  
  ## set condition
  N <- condition$N[cc]
  L <- condition$L[cc]
  K <- condition$K[cc]
  NN <- condition$NN[cc]
  VAR <- condition$VAR[cc]
  
  
  
  ## Result vectors
  Z_est <- C_est <- list()
  beta0_est <- beta1_est <- alpha_est <- gamma_est <- list()
  mu_est <- lambda_est <- B_est <- list()
  TIE_est <- DE_est <- rep(NA, reps)
  TIE_CI_est <- DE_CI_est <- matrix(NA, nrow = reps, ncol = 2)
  true_omegas <- array(NA, dim = c(N, L, reps))
  selected_list <- list()
  
  
  
  
  ## set true parameter values
  # class specific means
  true_mus <- true_mus_list[[cc]]
  # class specific covariance
  true_lambdas <- seq(1, 1 * VAR, length.out = L) * VAR
  true_B_temp <- seq(1, 1.2, length.out = K)
  true_B <- diag(true_B_temp / prod(true_B_temp)^(1 / K)) # normalize by the geometric mean
  true_Sigmas <- array(NA, dim = c(K, K, L))
  for (l in 1:L) {
    true_Sigmas[, , l] <- true_lambdas[l] * true_B
  }
  # regression coefficients gamma(LR), alpha (LR), beta (MLR)
  true_gamma <- 0
  true_alpha <- seq(-1, 1, length.out = L)
  true_beta1 <- seq(0, -L * 0.5, length.out = L)
  true_beta0 <- seq(0, 1, length.out = L)
  # calculate true classification probabilities
  pis <- matrix(NA, nrow = L, ncol = 2)
  for(l in 1:L){
    omega_temp <- rep(0, L)
    omega_temp[l] <- 1
    pis[l, 1] <- multinom_prob(omega_temp, true_beta0, true_beta1, 1)
    pis[l, 2] <- multinom_prob(omega_temp, true_beta0, true_beta1, 0)
  }
  
  
  
  
  ## start replication
  for (r in 1:reps) {
    
    set.seed(10 * as.numeric(args[1]) + r)
    
    # data generation
    data_temp <- sim_lcma(N = N, 
                          true_mus = true_mus, 
                          true_lambdas = true_lambdas, true_B = true_B, 
                          true_beta0 = true_beta0, true_beta1 = true_beta1,
                          true_alpha = true_alpha, true_gamma = true_gamma)
    true_omega <- data_temp$true_omega # true latent classes
    data <- add_noise(data_temp, NN)
    # plot(data$M[, c(1, 2)], col = apply(true_omega, 1, function(x){which(x == 1)}), main = paste("Condition", cc))
    # legend("topright", legend = paste0("Class", 1:4), col = 1:4, pch = 1)
    
    
    
    ## Variable Selection
    # var_select_TIE <- va1r_select(data, B = 1000, cr = "TIE")
    # var_select_DE <- var_select(data, B = 1000, cr = "DE")
    var_select_em_TIE <- var_select(data, B = 1000, cr = "em_TIE")
    
    # var_select_em_DE <- var_select(data, B = 1000, cr = "em_DE")
    selected <- var_select_em_TIE$selected[[length(var_select_em_TIE$selected)]]
    
    
    # Final model result
    result <- var_select_em_TIE$fit_final
    # relabel classes to prevent label switching by maximizing classification accuracy
    if (ncol(result$est$C) == L){
      C_temp <- result$est$C
      labels <- permn(L)
      class_accuracy <- sapply(labels, function(x) {
        mean(rowMeans(C_temp[, x] == true_omega) == 1)
        }
        )
      result$est$C <- C_temp[, labels[[which.max(class_accuracy)]]]
    }
    
    
    
    
    ## Bootstrap CI
    # boot_sample <- bootstrap(result, B = 1000)
    boot_sample <- emp_bootstrap(result, B = 1000)
    TIE_CI <- quantile(boot_sample[, "TIE"], c(0.025, 0.975), na.rm = TRUE)
    DE_CI <- quantile(boot_sample[, "DE"], c(0.025, 0.975), na.rm = TRUE)
    
    
    
    
    ## save results
    Z_est[[r]] <- result$est$Z
    C_est[[r]] <- result$est$C
    beta0_est[[r]] <- result$est$beta0
    beta1_est[[r]] <- result$est$beta1
    alpha_est[[r]] <- result$est$alpha
    gamma_est[[r]] <- result$est$gamma
    mu_est[[r]] <- result$est$mu
    lambda_est[[r]] <- result$est$lambda
    B_est[[r]] <- result$est$B
    
    
    true_omegas[, , r] <- data$true_omega
    TIE_est[r] <- result$est$TIE
    DE_est[r] <- result$est$DE
    TIE_CI_est[r, ] <- TIE_CI
    DE_CI_est[r, ] <- DE_CI
    selected_list[[r]] <- selected
    
    
    cat("Con", cc, "Rep", r, "\n")
    print(Sys.time() - start)
    
    
  } # end of replications
  
  
  parm_est <- list(Z = Z_est,
                   C = C_est,
                   beta0 = beta0_est,
                   beta1 = beta1_est,
                   alpha = alpha_est,
                   gamma = gamma_est,
                   mu = mu_est,
                   lambda = lambda_est,
                   B = B_est,
                   TIE = TIE_est,
                   DE = DE_est,
                   TIE_CI = TIE_CI_est,
                   DE_CI = DE_CI_est,
                   selected_list = selected_list)
  
  
  true_parameter <- list(mus = true_mus, 
                         lambdas = true_lambdas,
                         B = true_B,
                         Sigmas = true_Sigmas,
                         gamma = true_gamma,
                         alpha = true_alpha,
                         beta0 = true_beta0,
                         beta1 = true_beta1,
                         omega = true_omegas,
                         pis = pis)
  
  
  saveRDS(parm_est, paste0('parm_est_con_', cc, "_batch", args[1], '.rds'))
  saveRDS(true_parameter, paste0('true_parameter_con_', cc, "_batch", args[1], '.rds'))
  saveRDS(selected_list, paste0('selected_list_con_', cc, "_batch", args[1], '.rds'))
  
  
  print(Sys.time() - start)
  
  
}



