# =============================================================================
# 03_gibbs_sampler.R
# Bayesian Statistics and MCMC - Exercise 1 (Question D)
#
# Custom Gibbs sampler for Bayesian linear regression with conjugate priors:
#
#   Model:    Y = X*beta + epsilon,  epsilon ~ N(0, sigma^2 * I_n)
#   Prior:    beta | sigma^2 ~ N_p(0, sigma^2 * M^{-1})
#             sigma^2 ~ IG(a, b)
#
# With M^{-1} = 10000 * I_p (diffuse prior), the posterior means
# are nearly identical to OLS -- the data dominate entirely.
#
# Full conditionals (derived analytically in Exercise 1B):
#   beta    | sigma^2, y, X  ~  N_p(beta*, sigma^2 * (M + X'X)^{-1})
#   sigma^2 | beta,    y, X  ~  IG((n+p)/2 + a,  b + 0.5*(||y-X*beta||^2 + beta'*M*beta))
#
# where beta* = (M + X'X)^{-1} (X'X) beta_hat  (posterior mean = shrinkage of OLS)
# =============================================================================

library(mvtnorm)  # for rmvnorm

source("01_data_simulation.R")

# =============================================================================
# Gibbs Sampler Function
# =============================================================================

gibbs_bayes_lm <- function(X, y,
                            M_diag  = 1 / 10000,  # diagonal value of M (not M^{-1})
                            a       = 0.0001,      # IG prior shape
                            b       = 0.0001,      # IG prior rate
                            n_iter  = 20000,
                            burn_in = 5000,
                            seed    = 123) {
  set.seed(seed)
  
  n <- nrow(X)   # number of observations (50)
  p <- ncol(X)   # number of parameters  (16: intercept + 15 predictors)
  
  # --- Prior precision matrix M = (1/10000) * I_p ---
  # M^{-1} = 10000 * I_p  =>  prior variance of each beta_j is 10000 * sigma^2
  # This is a very diffuse prior; data will dominate
  M <- diag(M_diag, p)
  
  # --- Precompute constants (computed once, outside the loop for efficiency) ---
  XtX      <- crossprod(X)              # X'X  (p x p)
  Xty      <- crossprod(X, y)           # X'y  (p x 1)
  V_beta   <- solve(M + XtX)            # (M + X'X)^{-1}: posterior covariance kernel
  beta_star <- as.numeric(V_beta %*% Xty)  # posterior mean of beta
  # Note: beta_star = (M + X'X)^{-1}(X'X) beta_hat
  # With M very small, beta_star ≈ beta_hat (OLS)
  
  # --- Storage ---
  beta_samples   <- matrix(0, nrow = n_iter, ncol = p)
  sigma2_samples <- numeric(n_iter)
  
  # --- Initial values ---
  beta_cur   <- rep(0, p)
  sigma2_cur <- 1.0
  
  # --- Gibbs sampling loop ---
  for (t in 1:n_iter) {
    
    # Step 1: Sample beta | sigma^2, y, X
    # beta | sigma^2, y, X ~ N_p(beta*, sigma^2 * V_beta)
    Sigma_beta <- sigma2_cur * V_beta
    beta_cur   <- as.numeric(rmvnorm(1, mean = beta_star, sigma = Sigma_beta))
    beta_samples[t, ] <- beta_cur
    
    # Step 2: Sample sigma^2 | beta, y, X
    # sigma^2 | beta, y, X ~ IG((n+p)/2 + a,  b + 0.5*(RSS + beta'*M*beta))
    resid         <- as.numeric(y - X %*% beta_cur)
    RSS           <- sum(resid^2)                    # ||y - X*beta||^2
    beta_TMbeta   <- M_diag * sum(beta_cur^2)        # beta'*M*beta (simplified: M = M_diag * I)
    shape_sig     <- (n + p) / 2 + a
    rate_sig      <- b + 0.5 * (RSS + beta_TMbeta)
    sigma2_cur    <- 1 / rgamma(1, shape = shape_sig, rate = rate_sig)
    sigma2_samples[t] <- sigma2_cur
  }
  
  # --- Discard burn-in ---
  keep <- (burn_in + 1):n_iter
  
  return(list(
    beta        = beta_samples[keep, ],
    sigma2      = sigma2_samples[keep],
    beta_star   = beta_star,    # analytical posterior mean (for verification)
    V_beta      = V_beta,       # posterior covariance kernel
    n_kept      = length(keep)
  ))
}

# =============================================================================
# Run Gibbs Sampler
# =============================================================================

X_full <- cbind(1, X)  # add intercept column -> 50 x 16

gibbs_out <- gibbs_bayes_lm(
  X       = X_full,
  y       = y,
  M_diag  = 1 / 10000,
  a       = 0.0001,
  b       = 0.0001,
  n_iter  = 20000,
  burn_in = 5000
)

# =============================================================================
# Results
# =============================================================================

beta_post_mean <- colMeans(gibbs_out$beta)
beta_post_sd   <- apply(gibbs_out$beta, 2, sd)
beta_ci        <- apply(gibbs_out$beta, 2, quantile, probs = c(0.025, 0.975))
sigma2_mean    <- mean(gibbs_out$sigma2)

param_names <- c("beta0", paste0("beta", 1:15))

results_gibbs <- data.frame(
  True       = true_beta,
  OLS        = round(ols_coef, 4),
  Bayes_Mean = round(beta_post_mean, 4),
  Bayes_SD   = round(beta_post_sd, 4),
  CI_2.5     = round(beta_ci[1, ], 4),
  CI_97.5    = round(beta_ci[2, ], 4),
  row.names  = param_names
)

cat("=== Gibbs Sampler Results ===\n")
print(results_gibbs)
cat("\nPosterior mean of sigma^2 (Gibbs):", round(sigma2_mean, 4), "\n")
cat("OLS estimate of sigma^2:           ", round(ols_sigma2, 4), "\n")
cat("True sigma^2:                      ", true_sigma2, "\n")

# --- Variable selection: 95% CrI excludes zero ---
excludes_zero <- beta_ci[1, ] > 0 | beta_ci[2, ] < 0
cat("\nVariables with 95% CrI excluding zero:\n")
cat(param_names[excludes_zero], "\n")

# =============================================================================
# Diagnostic Plots
# =============================================================================

par(mfrow = c(3, 2))

# Trace plots
plot(gibbs_out$beta[, 1], type = "l", col = "steelblue",
     main = "Trace plot: beta0", xlab = "Iteration", ylab = expression(beta[0]))
plot(gibbs_out$beta[, 2], type = "l", col = "steelblue",
     main = "Trace plot: beta1", xlab = "Iteration", ylab = expression(beta[1]))
plot(gibbs_out$beta[, 6], type = "l", col = "steelblue",
     main = "Trace plot: beta5", xlab = "Iteration", ylab = expression(beta[5]))
plot(gibbs_out$beta[, 8], type = "l", col = "steelblue",
     main = "Trace plot: beta7", xlab = "Iteration", ylab = expression(beta[7]))
plot(gibbs_out$sigma2,    type = "l", col = "steelblue",
     main = "Trace plot: sigma2", xlab = "Iteration", ylab = expression(sigma^2))

par(mfrow = c(1, 1))

# ACF plots
par(mfrow = c(3, 2))
acf(gibbs_out$beta[, 1],  main = "ACF: beta0",  lag.max = 50)
acf(gibbs_out$beta[, 2],  main = "ACF: beta1",  lag.max = 50)
acf(gibbs_out$beta[, 6],  main = "ACF: beta5",  lag.max = 50)
acf(gibbs_out$beta[, 8],  main = "ACF: beta7",  lag.max = 50)
acf(gibbs_out$sigma2,     main = "ACF: sigma2", lag.max = 50)
par(mfrow = c(1, 1))

# Posterior boxplots with true values
boxplot(gibbs_out$beta,
        names  = param_names,
        col    = "orange",
        main   = "Posterior Distributions of Regression Coefficients (Gibbs)",
        xlab   = "Coefficient",
        ylab   = "Value",
        las    = 2)
abline(h = 0, lty = 2, col = "gray50")
points(1:16, true_beta, pch = 18, col = "red", cex = 1.5)
legend("topright", legend = "True value", pch = 18, col = "red")
