# =============================================================================
# 04_metropolis_hastings.R
# Bayesian Statistics and MCMC - Exercise 4
#
# Metropolis-Hastings sampling for the posterior of lambda under the
# Jeffreys prior, using data from Exercise 2 (T=35 errors in n=10 days).
#
# Model:   X_i ~ Poisson(lambda),  i = 1, ..., 10
#          T = sum(X_i) = 35  ~  Poisson(10 * lambda)
#
# Jeffreys prior:  p(lambda) ∝ lambda^{-1/2}   (improper, but posterior is proper)
#
# Posterior (analytical):  lambda | T ~ Gamma(T + 0.5, n) = Gamma(35.5, 10)
#   - Derived by: posterior ∝ likelihood x prior
#                           ∝ lambda^35 * exp(-10*lambda) * lambda^{-0.5}
#                           = lambda^{34.5} * exp(-10*lambda)
#                           => Gamma(35.5, 10)
#
# The MH algorithm is used to verify this analytically derived result.
# Proposal: Normal random walk  lambda* ~ N(lambda_current, sigma^2 = 8)
# Acceptance rate ~25% as required by the assignment.
# =============================================================================

T_obs <- 35   # total observed errors
n_obs <- 10   # number of days

# =============================================================================
# Log-posterior (unnormalized)
# We work on the log scale for numerical stability.
# The normalizing constant cancels in the MH ratio, so we only need the kernel.
# =============================================================================

log_posterior <- function(lambda) {
  if (lambda <= 0) return(-Inf)  # lambda must be positive; auto-reject if proposed <= 0
  (T_obs - 0.5) * log(lambda) - n_obs * lambda
  # = 34.5 * log(lambda) - 10 * lambda   (log of Gamma(35.5, 10) kernel)
}

# =============================================================================
# Metropolis-Hastings Function
# =============================================================================

metropolis_hastings <- function(n_iter, sigma2, seed = 123) {
  set.seed(seed)
  
  samples <- numeric(n_iter)
  accept  <- 0
  
  # Start from MLE = T/n = 3.5 (close to posterior mean 3.55 => fast convergence)
  current <- T_obs / n_obs
  
  for (i in 1:n_iter) {
    
    # --- Propose new value: Normal random walk ---
    candidate <- rnorm(1, mean = current, sd = sqrt(sigma2))
    
    # --- Log acceptance ratio ---
    # Because the Normal proposal is symmetric: q(lambda*|lambda) = q(lambda|lambda*)
    # the proposal terms cancel, leaving only the posterior ratio.
    # The normalizing constant of the posterior also cancels.
    log_alpha <- log_posterior(candidate) - log_posterior(current)
    
    # --- Accept / Reject ---
    # Equivalent to: accept if U < alpha, where U ~ Uniform(0,1)
    # Working on log scale avoids numerical underflow for very small alpha
    if (log(runif(1)) < log_alpha) {
      current <- candidate
      accept  <- accept + 1
    }
    
    # Store current state (whether accepted or not)
    samples[i] <- current
  }
  
  return(list(
    samples     = samples,
    accept_rate = accept / n_iter
  ))
}

# =============================================================================
# Tune sigma^2: find value giving ~25% acceptance rate
# =============================================================================

cat("=== Tuning sigma^2 for ~25% acceptance rate ===\n")
sigma2_grid <- c(0.1, 0.5, 1, 2, 3, 5, 8, 10, 15)

for (s2 in sigma2_grid) {
  test <- metropolis_hastings(20000, s2)
  cat(sprintf("sigma^2 = %5.1f  ->  acceptance rate = %.3f\n",
              s2, test$accept_rate))
}

# sigma^2 = 8 gives ~25% acceptance rate
sigma2_opt <- 8

# =============================================================================
# Run MH with optimal sigma^2
# =============================================================================

total_iter  <- 21000
burnin      <- 1000

result <- metropolis_hastings(total_iter, sigma2_opt)
cat(sprintf("\nFinal acceptance rate (sigma^2 = %g): %.3f\n",
            sigma2_opt, result$accept_rate))

# Discard burn-in
lambda_samples <- result$samples[(burnin + 1):total_iter]
cat("Kept", length(lambda_samples), "samples after burn-in of", burnin, "\n")

# =============================================================================
# Results: compare MH samples with analytical Gamma(35.5, 10)
# =============================================================================

# Analytical posterior moments
alpha_post <- T_obs + 0.5   # = 35.5
beta_post  <- n_obs          # = 10
analytical_mean <- alpha_post / beta_post           # = 3.55
analytical_sd   <- sqrt(alpha_post / beta_post^2)  # ≈ 0.596
analytical_ci   <- qgamma(c(0.025, 0.975), shape = alpha_post, rate = beta_post)

# MH sample moments
mh_mean <- mean(lambda_samples)
mh_sd   <- sd(lambda_samples)
mh_ci   <- quantile(lambda_samples, probs = c(0.025, 0.975))

cat("\n=== Comparison: MH vs Analytical Gamma(35.5, 10) ===\n")
cat(sprintf("%-20s %10s %10s\n", "", "MH", "Analytical"))
cat(sprintf("%-20s %10.4f %10.4f\n", "Mean",   mh_mean, analytical_mean))
cat(sprintf("%-20s %10.4f %10.4f\n", "SD",     mh_sd,   analytical_sd))
cat(sprintf("%-20s %10.4f %10.4f\n", "2.5%",  mh_ci[1], analytical_ci[1]))
cat(sprintf("%-20s %10.4f %10.4f\n", "97.5%", mh_ci[2], analytical_ci[2]))

cat("\n=== 95% Credible Interval (MH) ===\n")
cat(sprintf("(%.3f, %.3f)\n", mh_ci[1], mh_ci[2]))
cat("Interpretation: Given T=35 errors in 10 days and the Jeffreys prior,\n")
cat("there is 95% probability that the true error rate lambda lies between",
    round(mh_ci[1], 3), "and", round(mh_ci[2], 3), "\n")
cat("Note: H0 value (lambda=2) is OUTSIDE this interval, consistent with\n")
cat("the strong evidence against H0 found in Exercise 2 (BF10 = 15.45).\n")

# =============================================================================
# Diagnostic Plots
# =============================================================================

par(mfrow = c(2, 2))

# 1. Trace plot
plot(lambda_samples, type = "l", col = "steelblue",
     main = "Trace Plot of lambda (MH algorithm)",
     xlab = "Iteration (after burn-in)", ylab = expression(lambda))
abline(h = analytical_mean, col = "red", lty = 2, lwd = 2)
legend("topright", legend = "Analytical mean", col = "red", lty = 2)

# 2. Ergodic mean plot
ergodic_mean <- cumsum(lambda_samples) / seq_along(lambda_samples)
plot(ergodic_mean, type = "l", col = "steelblue",
     main = "Ergodic Mean Plot of lambda",
     xlab = "Iteration (after burn-in)", ylab = expression(lambda))
abline(h = analytical_mean, col = "red", lty = 2, lwd = 2)
legend("topright", legend = sprintf("Analytical mean = %.2f", analytical_mean),
       col = "red", lty = 2)

# 3. ACF plot
acf(lambda_samples, lag.max = 50,
    main = "Autocorrelation of lambda (MH algorithm)")

# 4. Posterior histogram vs analytical density
hist(lambda_samples, breaks = 40, probability = TRUE,
     col  = "lightblue", border = "white",
     main = "Posterior Distribution of lambda (MH)",
     xlab = expression(lambda))
curve(dgamma(x, shape = alpha_post, rate = beta_post),
      add = TRUE, col = "red", lwd = 2)
abline(v = mh_ci, col = "darkgreen", lty = 2, lwd = 1.5)
legend("topright",
       legend = c("Analytical Gamma(35.5, 10)", "95% CI"),
       col    = c("red", "darkgreen"),
       lty    = c(1, 2), lwd = c(2, 1.5))

par(mfrow = c(1, 1))
