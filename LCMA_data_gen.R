## function for the data generation of the latent class mediation analysis


sim_lcma <- function(N, # sample size
                     true_mus, true_lambdas, true_B, 
                     true_beta0, true_beta1, 
                     true_alpha, true_gamma) {
  
  L <- ncol(true_mus)
  K <- nrow(true_mus)
  
  # true covariance matrix
  true_Sigmas <- array(NA, dim = c(K, K, L))
  for (l in 1:L) {
    true_Sigmas[, , l] <- true_lambdas[l] * true_B
  }
  # Group membership
  G <- rbinom(N, 1, 0.5)
  # true classification probabilities
  pis <- matrix(NA, nrow = L, ncol = 2)
  for(l in 1:L){
    omega_temp <- rep(0, L)
    omega_temp[l] <- 1
    pis[l, 1] <- multinom_prob(omega_temp, true_beta0, true_beta1, 1)
    pis[l, 2] <- multinom_prob(omega_temp, true_beta0, true_beta1, 0)
  }
  # true class membership
  true_omega <- matrix(0, nrow = N, ncol = L)
  for (i in 1:N){
    true_omega[i, sample(x = 1:L, 1, prob = pis[, 2 - G[i]])] <- 1
  }
  # Response data Y
  Y <- rep(NA, N)
  for (i in 1:N) {
    Y[i] <- rbinom(1, 1, logistic_prob(true_gamma, true_alpha, true_omega[i, ], G[i]))
  }
  # Process feature M
  M <- matrix(NA, nrow = N, ncol = K)
  for (l in 1:L) {
    M[true_omega[, l] == 1, ] <- MASS::mvrnorm(sum(true_omega[, l] == 1),
                                               mu = true_mus[, l],
                                               Sigma = true_Sigmas[, , l])
  }
  
  
  inputs = list(N = N,
                true_mus = true_mus,
                true_lambdas = true_lambdas,
                true_B = true_B,
                true_beta0 = true_beta0,
                true_beta1 = true_beta1,
                true_alpha = true_alpha,
                true_gamma = true_gamma)
  
  
  data <- list(Y = Y,
               M = M,
               G = G,
               true_omega = true_omega,
               pis = pis,
               inputs = inputs)
  
  
  data
}



