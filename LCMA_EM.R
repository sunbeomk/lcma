## function for latent class mediation analysis using EM algorithm
library(MASS)
library(emdbook)




## EM algorithm for latent class mediation analysis
lcma <- function(Y,    # binary response vector
                 M,    # continuous process feature
                 G,    # binary group membership vector
                 L,    # number of latent classes
                 Z0,   # initial guess for class membership
                 epsilon = 1e-4, # tolerance for EM algorithm
                 m = 100L        # maximum number of iteration
                 ) {
  
  M <- as.matrix(M)
  
  K <- ncol(M)
  N <- length(Y)
  
  
  ## setting parameter vectors
  Z <- array(NA, dim = c(N, L, m))
  C <- array(NA, dim = c(N, L, m))
  beta0 <- beta1 <- alpha <- lambda <- matrix(NA, nrow = L, ncol = m)
  gamma <- rep(NA, m)
  mu <- array(NA, dim = c(K, L, m))
  B <- array(NA, dim = c(K, K, m))
  class_accuracy <- rep(NA, m)
  class_no <- matrix(NA, nrow = N, ncol = m)
  
  ## initial values
  Z[, , 1] <- C[, , 1] <- Z0
  if (K > 1) {
    B[, , 1] <- diag(diag(cov(M))) / det(cov(M))^{1 / K}
  } else {
    B[, , 1] <- var(M) / det(var(M))^{1 / K}
  }
  beta0_iter <- beta1_iter <- rep(0, L)
  gamma_iter <- 0
  alpha_iter <- rep(0, L)
  
  ## EM algorithm
  iter <- 1L   # iteration
  loglik_vec <- rep(NA, m)
  start <- Sys.time()
  while (TRUE) {
    # M-step ----
    ## update beta ----
    betas_temp <- optim(par = c(beta0_iter[-1], beta1_iter[-1]), 
                        fn = neg_multinom_loglik, 
                        Z = Z[, , iter], 
                        G = G, 
                        method = "BFGS", 
                        control = list(ndeps = rep(1e-8, 2 * (L - 1)),
                                       maxit = 1000))$par
    beta0[, iter] <- beta0_iter <- c(0, betas_temp[1:(L - 1)])
    beta1[, iter] <- beta1_iter <- c(0, betas_temp[L:(2 * (L - 1))])
    ## update gamma, alpha ----
    gamma_alpha_temp <- optim(par = c(gamma_iter, alpha_iter), 
                              fn = neg_logistic_loglik, 
                              Z = Z[, , iter], 
                              G = G, 
                              Y = Y, 
                              method = "BFGS",
                              control = list(ndeps = rep(1e-8, 1 + L),
                                             maxit = 1000))$par
    gamma[iter] <- gamma_iter <- gamma_alpha_temp[1]
    alpha[, iter] <- alpha_iter <- gamma_alpha_temp[-1]
    ## update mu ----
    n_ls <- colSums(Z[, , iter])
    mu[, , iter] <- t(apply(t(M) %*% Z[, , iter], 1, function(x) {x / n_ls}))
    ## update Sigma ----
    W_l <- array(NA, dim = c(K, K, L))    # scattering matrix
    B_temp <- array(NA, dim = c(K, K, L))
    for (l in 1:L) {
      # W_l[, , l] <- cov(M[which(C[, l, iter] == 1), ]) * n_ls[l]
      W_l_temp <- array(NA, dim = c(K, K, N))
      weighted_means <- colSums(apply(as.matrix(M), 2, function(x) {x * Z[, l, iter]})) / n_ls[l]
      for (i in 1:N) {
        W_l_temp[, , i] <- Z[i, l, iter] * ((as.matrix(M)[i, ] - weighted_means) %*% t(as.matrix(M)[i, ] - weighted_means))
      }
      W_l[, , l] <- apply(W_l_temp, 2, rowSums)
    }
    
    
    ## update lambda and B iteratively ----
    for (mm in 1:50){
      for (l in 1:L) {
        lambda[l, iter] <- sum(diag(W_l[, , l] %*% solve(B[, , iter]))) / (K * n_ls[l])
        B_temp[, , l] <- (1 / lambda[l, iter]) * W_l[, , l]
      }
      W_l_sum <- apply(B_temp, 1, rowSums)
      if (K > 1) {
        B[, , iter] <- diag(diag(W_l_sum) / det(W_l_sum)^(1 / K))
      } else {
        B[, , iter] <- W_l_sum / det(as.matrix(W_l_sum))^(1 / K)
      }
    }
    
    
    ## model log-likelihood
    loglik_vec[iter] <- model_loglik(mu[ , ,iter], lambda[, iter], B[, , iter], 
                                     gamma[iter], alpha[, iter], 
                                     beta0[, iter], beta1[, iter], 
                                     Y, M, G)
    # class_accuracy[iter] <- class_acc(C[,,iter], true_omega)
    class_no[, iter] <- apply(C[,,iter], 1, which.max)
    
    
    ## termination
    if (iter == m) {break}
    if (iter > 1) {
      if ((loglik_vec[iter] - loglik_vec[iter - 1]) < epsilon) {
        iter <- iter - 1L
        break}
    }
    iter <- iter + 1L
    # set initial B for the next iteration
    B[, , iter] <- B[, , iter - 1]
    
    
    
    
    # E-step ----
    # update Z
    # for (i in 1:N) {
    #   P_OYM <- rep(NA, L)
    #   for (l in 1:L) {
    #     # P_omega
    #     P_omega <- exp(beta0[l, iter - 1] + beta1[l, iter - 1] * G[i]) / sum(exp(beta0[, iter - 1] + beta1[, iter - 1] * G[i]))
    #     # P_Y
    #     exp_Y <- exp(gamma[iter - 1] * G[i] + alpha[l, iter - 1])
    #     P_correct <- exp_Y / (1 + exp_Y)
    #     P_Y <- Y[i] * P_correct + (1 - Y[i]) * (1 - P_correct)
    #     # P_M
    #     P_M <- dmvnorm(M[i, ], mu[, l, iter - 1], lambda[l, iter - 1] * B[, , iter - 1])
    #     
    #     P_OYM[l] <- P_omega * P_Y * P_M
    #   }
    #   
    #   Z[i, , iter] <- P_OYM / sum(P_OYM)
    # }
    for (g in c(0, 1)) {
      for (y in c(0, 1)) {
        P_OYM <- matrix(NA, nrow = sum(G == g & Y == y), ncol = L)
        for (l in 1:L) {
          # P_omega
          P_omega <- exp(beta0[l, iter - 1] + beta1[l, iter - 1] * g) / sum(exp(beta0[, iter - 1] + beta1[, iter - 1] * g))
          # P_Y
          exp_Y <- exp(gamma[iter - 1] * g + alpha[l, iter - 1])
          P_correct <- exp_Y / (1 + exp_Y)
          P_Y <- y * P_correct + (1 - y) * (1 - P_correct)
          # P_M
          P_M <- apply(as.matrix(M[G == g & Y == y, ]), 1, dmvnorm, mu = mu[, l, iter - 1], sigma = lambda[l, iter - 1] * B[, , iter - 1])
          # P_OYM
          P_OYM[, l] <- P_omega * P_Y * P_M
        }
        
        Z[G == g & Y == y, , iter] <- apply(P_OYM, 1, function(x){x / sum(x)})
      }
    }
    
    
    # update C
    C[, , iter] <- classify_Z(Z[, , iter])
    
    
    # cat("Iteration", iter, "")
    # print(Sys.time() - start)
    ## End of EM algorithm
    
    
  }
  
  
  ## Total Indirect Effect, Direct Effect
  tie <- TIE(beta0[, iter], beta1[, iter], gamma[iter], alpha[, iter])
  de <- DE(beta0[, iter], beta1[, iter], gamma[iter], alpha[, iter])
  
  
  ## save results
  est <- list(Z = Z[, , iter],
              C = C[, , iter],
              beta0 = beta0[, iter],
              beta1 = beta1[, iter],
              alpha = alpha[, iter],
              gamma = gamma[iter],
              mu = mu[, , iter],
              lambda = lambda[, iter],
              B = B[, , iter],
              TIE = tie,
              DE = de)
  
  
  parm <- list(Z = Z[, , 1:iter],
               C = C[, , 1:iter],
               beta0 = beta0[, 1:iter],
               beta1 = beta1[, 1:iter],
               alpha = alpha[, 1:iter],
               gamma = gamma[1:iter],
               mu = mu[, , 1:iter],
               lambda = lambda[, 1:iter],
               B = B[, , 1:iter],
               class_accuracy = class_accuracy[1:iter])
  
  
  inputs = list(Y = Y,
                M = M,
                G = G,
                L = L)
  
  
  result <- list(est = est,
                 inputs = inputs,
                 parm = parm,
                 iter = iter,
                 loglik_vec = loglik_vec[1:iter]
                 )
  
  
  result
}




# result <- lcma(Y, M, G, L, Z0)
# result <- lcma(Y = data_temp$Y, 
#                M = data_temp$M, 
#                G = data_temp$G, 
#                L = L, 
#                Z0 = Z0)




# # classification accuracy
# mean(rowMeans(classify_Z(result$est$Z) == true_omega) == 1)
# model implied classification probabilities
# model_implied_pis <- matrix(NA, nrow = L, ncol = 2)
# for(l in 1:L){
#   omega_temp <- rep(0, L)
#   omega_temp[l] <- 1
#   model_implied_pis[l, 1] <- multinom_prob(omega_temp,
#                                            beta0[, iter],
#                                            beta1[, iter], 1)
#   model_implied_pis[l, 2] <- multinom_prob(omega_temp,
#                                            beta0[, iter],
#                                            beta1[, iter], 0)
# }
# model_implied_pis # estimated probabilities
# pis # true probabilities
# (model_implied_pis - pis) / pis # relative bias




# class_accuracy <- rep(NA, iter)
# for (i in 1:iter) {
#   class_accuracy[i] <- mean(rowMeans(C[,,i] == true_omega) == 1)
# }
# plot(class_accuracy, type = "l")
# 
# 
# plot(loglik_vec[1:iter], type = "l")
# 
# 
# matplot(t(beta0[, 1:iter]), type = "l", lty = 1, ylim = range(true_beta0))
# abline(h = true_beta0, lty = 2, col = 1:L)
# 
# 
# matplot(t(beta1[, 1:iter]), type = "l", lty = 1, ylim = range(true_beta1))
# abline(h = true_beta1, lty = 2, col = 1:L)
# 
# 
# plot(gamma, type = "l", ylim = c(true_gamma - 0.1, true_gamma + 0.1))
# 
# 
# matplot(t(alpha[, 1:iter]), type = "l", lty = 1, ylim = range(true_alpha))
# abline(h = true_alpha, lty = 2, col = 1:L)
# 
# 
# mu_mse <- rep(NA, iter)
# for (i in 1:iter) {
#   mu_mse[i] <- mean((mu[, , i] - true_mus)^2)
# }
# plot(mu_mse, type = "l")
# 
# 
# matplot(t(lambda[, 1:iter]), type = "l", lty = 1, ylim = c(0, 5))
# abline(h = true_lambdas, lty = 2, col = 1:L)
# 
# 
# B_mse <- rep(NA, iter)
# for (i in 1:iter) {
#   B_mse[i] <- mean((diag(B[, , i]) - diag(true_B))^2)
# }
# plot(B_mse, type = "l")
# 
# 
# diag(B[,,iter])
# diag(true_B)
# 
# 
# true_class <- apply(true_omega, 1, which.max)
# for (i in 1:iter) {
#   class_no <- apply(C[, , i], 1, which.max)
#   print(i)
#   print(table(class_no, true_class))
# }



