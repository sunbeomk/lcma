library(parallel)


# multivariate normal density ----
dmvnorm <- function(x, mu, sigma) {
  k <- length(x)
  if (k == 1) {
    sigma <- as.matrix(sigma)
    mu <- t(as.matrix(mu))
  }
  (2 * pi)^{-k / 2} * 
    det(sigma)^{-1 / 2} * 
    exp((-1 / 2) * t(x - mu) %*% solve(sigma) %*% (x - mu))
}


# multinomial probability ----
multinom_prob <- function(omega_i, beta0, beta1, g_i) {
  exp((t(beta0) %*% omega_i + (t(beta1) %*% omega_i) * g_i)) / sum(exp(beta0 + beta1 * g_i))
}


# logistic probability ----
logistic_prob <- function(gamma, alpha, omega_i, g_i) {
  expo <- exp(gamma * g_i + t(alpha) %*% omega_i)
  expo / (1 + expo)
}


# -1 * multinomial log-likelihood ----
neg_multinom_loglik <- function(parm, Z, G) {
  L <- ncol(Z)
  N <- length(G)
  
  beta0 <- c(0, parm[1:(L - 1)])
  beta1 <- c(0, parm[L:(2 * (L - 1))])
  

  log_prob <- matrix(NA, nrow = N, ncol = L)
  for (l in 1:L) {
    for (g in c(0, 1)) {
      log_prob[which(G == g), l] <- Z[G == g, l] * log(exp(beta0[l] + beta1[l] * g) / sum(exp(beta0 + beta1 * g)))
    }
  }
  -1 * sum(log_prob)
}


# -1 * logistic log-likelihood ----
neg_logistic_loglik <- function(parm, Z, G, Y) {
  gamma <- parm[1]
  alpha <- parm[-1]
  L <- ncol(Z)
  
  log_prob <- matrix(NA, nrow = N, ncol = L)
  for (l in 1:L) {
    for (g in c(0, 1)) {
      correct_prob <- exp(gamma * g + alpha[l]) / (1 + exp(gamma * g + alpha[l]))
      log_prob[which(G == g), l] <- Z[G == g, l] * (Y[G == g] * log(correct_prob) + (1 - Y[G == g]) * log(1 - correct_prob))
      
    }
  }
  
  -1 * sum(log_prob)
}


# model log-likelihood ----
model_loglik <- function(mu, lambda, B, gamma, alpha, beta0, beta1, Y, M, G) {
  
  if (length(M) == length(Y)) {
    B <- as.matrix(B)
    mu <- t(as.matrix(mu))
    M <- as.matrix(M)
    
    L <- length(mu)
    K <- 1
  } else {
    K <- ncol(M)
    L <- ncol(mu)
  }
  N <- length(Y)
  
  likelihood_temp <- matrix(NA, nrow = N, ncol = L)
  for (i in 1:N) {
    for(l in 1:L) {
      multinom <- exp(beta0[l] + beta1[l] * G[i]) / sum(exp(beta0 + beta1 * G[i]))
      
      prob_Y <- exp(gamma * G[i] + alpha[l]) / (1 + exp(gamma * G[i] + alpha[l]))
      logistic <- Y[i] * prob_Y + (1 - Y[i]) * prob_Y
      
      multinorm <- dmvnorm(M[i, ], mu[, l], lambda[l] * B)
      
      # mclust::dmvnorm(t(M[i, ]), mu[, l], lambda[l] * B, log = TRUE)
      likelihood_temp[i, l] <- prod(multinom, logistic, multinorm)
    }
  }

  sum(log(rowSums(likelihood_temp)))
}


# classification from Z ----
classify_Z <- function(Z) {
  L <- ncol(Z)
  N <- nrow(Z)
  
  Z_class <- matrix(0, nrow = N, ncol = L)
  
  for (i in 1:N) {
    Z_class[i, which.max(Z[i, ])] <- 1
  }
  
  Z_class
}


# classification accuracy ----
class_acc <- function(estimate, true) {
  mean(rowMeans(estimate == true) == 1)
}


# contingency table ----
contingency <- function(estimate, true) {
  class_no <- apply(estimate, 1, which.max)
  true_no <- apply(true, 1, which.max)
  
  table(class_no, true_no)
}


# the expectation of Y in group g, when the latent class mediator \Omega is held
# constant at teh value it would obtain for gorup g'
E_Y <- function(g, gp, beta0, beta1, gamma, alpha) {
  L <- length(beta0)
  
  P_Y <- exp(gamma * g + alpha) / (1 + exp(gamma * g + alpha))
  E_Y <- 1 * P_Y + 0 * (1 - P_Y)
  
  P_O <- exp(beta0 + beta1 * gp) / sum(exp(beta0 + beta1 * gp))
  
  sum(E_Y * P_O)
}


# Total Indirect Effect ----
TIE <- function(beta0, beta1, gamma, alpha) {
  E_Y(g = 1, gp = 1, beta0, beta1, gamma, alpha) - 
    E_Y(g = 1, gp = 0, beta0, beta1, gamma, alpha)
}


# Direct Effect ----
DE <- function(beta0, beta1, gamma, alpha) {
  E_Y(g = 1, gp = 0, beta0, beta1, gamma, alpha) - 
    E_Y(g = 0, gp = 0, beta0, beta1, gamma, alpha)
}


# Bootstrap with a single sample ----
boot <- function(x, result) {
  set.seed(x)
  TIE_value <- DE_value <- NA
  # input
  Y <- result$inputs$Y
  G <- result$inputs$G
  M <- result$inputs$M
  C <- result$est$C
  
  N <- nrow(M)
  K <- ncol(M)
  L <- result$inputs$L
  # bootstrap sample
  ind <- sample(1:N, N, replace = TRUE)
  # single indicator or multiple indicator
  model <- if (K > 1) {"VEI"} else {NULL}
  # initial guessing for Z
  mclust_fit <- mclust::Mclust(data = M[ind, ], 
                               G = L, 
                               verbose = FALSE,
                               modelNames = model)
  if (is.null(mclust_fit)) {return(c(NA, NA))}
  Z0_no <- mclust_fit$classification
  Z0_temp <- matrix(0, nrow = N, ncol = L)
  for (i in 1:N) {
    Z0_temp[i, Z0_no[i]] <- 1
  }
  Z0 <- Z0_temp
  
  # run LCMA
  tryCatch({  
    result_temp <- lcma(Y = Y[ind], M = M[ind, ], G = G[ind], L = L, Z0 = Z0)
    # save result
    TIE_value <- result_temp$est$TIE
    DE_value <- result_temp$est$DE
  }, error = function(e) NA)
  
  # # run LCMA
  # result_temp <- lcma(Y = Y[ind], M = M[ind, ], G = G[ind], L = L, Z0 = Z0)
  # # save result
  # TIE <- result_temp$est$TIE
  # DE <- result_temp$est$DE
  
  c(TIE_value, DE_value)
}


# Bootstrap Confidence Interval for Total Indirect Effect / Direct Effect ----
bootstrap <- function(result, B = 1000) {
  # run bootstrap
  numCores <- min(8, detectCores())
  boot_result <- mclapply(1:B, boot, mc.cores = numCores, result = result)
  # save results as a matrix
  result <- matrix(unlist(boot_result), ncol = 2, byrow = TRUE)
  colnames(result) <- c("TIE", "DE")
  result
}


# Add NN noisy variables to the data ----
add_noise <- function(data, NN) {
  ## NN: number of noisy features
  
  N <- data$inputs$N
  L <- length(data$inputs$true_lambdas)
  K <- ncol(data$M)
  
  M_temp <- MASS::mvrnorm(n = N,
                          mu = rep(-10, NN),
                          Sigma = diag(runif(NN, min(data$inputs$true_lambdas) * min(diag(data$inputs$true_B)), 
                                                 max(data$inputs$true_lambdas) * max(diag(data$inputs$true_B))))
  )
  
  # replace with noisy variables
  data$M[, (K - NN + 1):K] <- M_temp
  
  data
}


# Calculate Empirical TIE and DE ----
emp_TIE <- function(Y, M, G) {
  
  N <- length(Y)
  K <- ncol(M)
  
  model <- if (K > 1) {"VEI"} else {NULL}
  mclust_fit <- mclust::Mclust(data = M, G = 2:6, verbose = FALSE, modelNames = model)
  # return NAs if mclust fails to fit
  if (is.null(mclust_fit)) {
    return(c(NA, NA))
    }
  L <- mclust_fit$G
  Z <- mclust_fit$classification
  
  # return NAs if empty class exists
  if (length(table(Z)) < L) {
    return(c(NA, NA))
  }
  
  p_y <- p_omega <- matrix(NA, nrow = L, ncol = 2)
  colnames(p_y) <- colnames(p_omega) <- c(1, 0)
  for (l in 1:L) {
    p_y[l, 1] <- mean(Y[G == 1 & Z == l])
    p_y[l, 2] <- mean(Y[G == 0 & Z == l])
  }
  
  for (l in 1:L) {
    p_omega[l, 1] <- sum(Z[G == 1] == l) / sum(G == 1)
    p_omega[l, 2] <- sum(Z[G == 0] == l) / sum(G == 0)
  }
  
  emp_TIE_est <- sum(p_y[, 1] * p_omega[, 1]) - sum(p_y[, 1] * p_omega[, 2])
  emp_DE_est <- sum(p_y[, 1] * p_omega[, 2]) - sum(p_y[, 2] * p_omega[, 2])
  
  c(emp_TIE_est, emp_DE_est)
}


# Bootstrap with a single sample ----
emp_boot <- function(x, result) {
  set.seed(x)
  TIE <- DE <- NA
  # input
  Y <- result$inputs$Y
  G <- result$inputs$G
  M <- result$inputs$M
  C <- result$est$C
  
  N <- nrow(M)
  K <- ncol(M)
  L <- result$inputs$L
  # bootstrap sample
  ind <- sample(1:N, N, replace = TRUE)
  # single indicator or multiple indicator
  
  emp_TIE(Y = Y[ind], M = as.matrix(M[ind, ]), G = G[ind])
  
}




# Bootstrap Confidence Interval for Total Indirect Effect / Direct Effect ----
emp_bootstrap <- function(result, B = 1000) {
  # run bootstrap
  numCores <- min(8, detectCores())
  boot_result <- mclapply(1:B, emp_boot, mc.cores = numCores, result = result)
  # save results as a matrix
  result <- matrix(unlist(boot_result), ncol = 2, byrow = TRUE)
  colnames(result) <- c("TIE", "DE")
  result
}



