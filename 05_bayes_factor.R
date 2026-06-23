# =============================================================================
# 05_bayes_factor.R
# Bayesian Statistics and MCMC - Exercise 2
#
# Bayesian hypothesis test for the Poisson error rate lambda:
#   H0: lambda = 2   vs   H1: lambda != 2
#
# Data: T = sum(X_i) = 35 errors observed over n = 10 days
#       T ~ Poisson(n * lambda) = Poisson(10 * lambda)
#
# Prior under H1: lambda ~ Gamma(shape=2, rate=1)
#   - Mean = 2/1 = 2  (coincides with H0 value => weak/fair prior)
#   - Variance = 2/1^2 = 2  (large relative to mean => diffuse belief)
#
# Prior probabilities: P(H0) = P(H1) = 0.5  (equally probable a priori)
#
# Bayes Factor: BF10 = P(T=35 | H1) / P(T=35 | H0)
# Posterior Odds: PO10 = BF10 * P(H1)/P(H0) = BF10 * 1 = BF10
#
# Posterior predictive distribution of X_{11} (next observation) under H1:
#   X_{11} | T ~ NegBin(r=37, p=11/12)
# =============================================================================

T_obs <- 35   # total observed errors
n_obs <- 10   # number of days

# =============================================================================
# Part A: Bayes Factor and Posterior Odds
# =============================================================================

# --- Marginal likelihood under H0 ---
# lambda is known (= 2), so T ~ Poisson(10*2) = Poisson(20)
# P(T=35 | H0) = exp(-20) * 20^35 / 35!
# We use log scale to avoid numerical overflow
log_m0 <- dpois(T_obs, lambda = n_obs * 2, log = TRUE)

# --- Marginal likelihood under H1 ---
# lambda is unknown with prior lambda ~ Gamma(2, 1)
# P(T=35 | H1) = integral_0^inf P(T=35 | lambda) * p(lambda) d_lambda
#
# Analytical solution (Gamma-Poisson conjugacy):
# integral_0^inf lambda^T * exp(-n*lambda) * lambda^(alpha-1) * exp(-beta*lambda) d_lambda
# = integral_0^inf lambda^(T+alpha-1) * exp(-(n+beta)*lambda) d_lambda
# = Gamma(T+alpha) / (n+beta)^(T+alpha)
#
# With T=35, alpha=2, beta=1, n=10:
# P(T=35 | H1) = 10^35 * 36 / 11^37   (using Gamma(37) = 36!)
log_m1 <- T_obs * log(10) + log(36) - 37 * log(11)
# Equivalent to: log(10^35) + log(36) - log(11^37)

# --- Bayes Factor and Posterior Odds ---
BF10 <- exp(log_m1 - log_m0)
# Because P(H0) = P(H1) = 0.5, prior odds = 1, so PO10 = BF10
PO10 <- BF10

cat("=== Bayesian Hypothesis Test: H0: lambda=2 vs H1: lambda!=2 ===\n\n")
cat(sprintf("log P(T=35 | H0) = %.4f\n", log_m0))
cat(sprintf("log P(T=35 | H1) = %.4f\n", log_m1))
cat(sprintf("Bayes Factor BF10 = %.4f\n", BF10))
cat(sprintf("Posterior Odds PO10 = %.4f\n\n", PO10))

# --- Jeffreys scale interpretation ---
cat("Interpretation (Jeffreys scale):\n")
if (BF10 < 3) {
  cat("  BF10 < 3: Anecdotal evidence for H1\n")
} else if (BF10 < 10) {
  cat("  3 < BF10 < 10: Moderate evidence for H1\n")
} else if (BF10 < 30) {
  cat("  10 < BF10 < 30: Strong evidence for H1\n")  # <-- our case
} else {
  cat("  BF10 > 30: Very strong evidence for H1\n")
}
cat(sprintf("  BF10 = %.2f: H1 is %.2f times more supported by the data than H0\n",
            BF10, BF10))
cat("  Conclusion: Reject H0. The true error rate lambda differs significantly from 2.\n")
cat(sprintf("  (Empirical rate: T/n = %d/%d = %.1f, well above lambda0 = 2)\n\n",
            T_obs, n_obs, T_obs / n_obs))

# =============================================================================
# Part B: Posterior Predictive Distribution of X_{11}
# =============================================================================

# --- Posterior of lambda under H1 (conjugate update) ---
# Prior: lambda ~ Gamma(alpha=2, beta=1)
# Likelihood: T ~ Poisson(n*lambda)
# Posterior: lambda | T ~ Gamma(alpha + T, beta + n) = Gamma(37, 11)
alpha_post <- 2 + T_obs   # = 37
beta_post  <- 1 + n_obs   # = 11

cat("=== Posterior of lambda under H1 ===\n")
cat(sprintf("lambda | T ~ Gamma(%d, %d)\n", alpha_post, beta_post))
cat(sprintf("Posterior mean:     %.4f  (prior mean was 2.0, MLE is %.1f)\n",
            alpha_post / beta_post, T_obs / n_obs))
cat(sprintf("Posterior variance: %.4f\n\n", alpha_post / beta_post^2))

# --- Posterior Predictive: X_{11} | T ~ NegBin(r, p) ---
# P(X_{11}=x | T) = integral P(X_{11}=x | lambda) * p(lambda | T) d_lambda
# = NegBin(r = alpha_post, p = beta_post / (beta_post + 1))
#
# This is the Poisson-Gamma mixture result:
# Poisson(lambda) x Gamma(r, beta) => marginal is NegBin(r, beta/(beta+1))
r_nb <- alpha_post                        # = 37
p_nb <- beta_post / (beta_post + 1)      # = 11/12

# Moments of NegBin(r, p): E = r*(1-p)/p,  Var = r*(1-p)/p^2
E_bayes   <- r_nb * (1 - p_nb) / p_nb
Var_bayes <- r_nb * (1 - p_nb) / p_nb^2

# Decomposition via Law of Total Variance:
# Var[X_{11}|T] = E[lambda|T] + Var[lambda|T]
E_lambda_post   <- alpha_post / beta_post       # = 37/11 ≈ 3.364
Var_lambda_post <- alpha_post / beta_post^2     # = 37/121 ≈ 0.306

# Classical (frequentist) predictive: X_{11} | lambda_hat ~ Poisson(lambda_hat)
lambda_mle  <- T_obs / n_obs   # = 3.5
E_class     <- lambda_mle
Var_class   <- lambda_mle

cat("=== Posterior Predictive Distribution of X_{11} ===\n")
cat(sprintf("Bayesian: X_{11} | T ~ NegBin(r=%d, p=%.4f)\n", r_nb, p_nb))
cat(sprintf("Classical: X_{11} | lambda_hat ~ Poisson(%.1f)\n\n", lambda_mle))

cat(sprintf("%-25s %10s %10s\n", "", "Bayesian", "Classical"))
cat(sprintf("%-25s %10.4f %10.4f\n", "Mean",     E_bayes,   E_class))
cat(sprintf("%-25s %10.4f %10.4f\n", "Variance", Var_bayes, Var_class))
cat("\nNote: Bayesian variance > Classical variance because:\n")
cat("  Var[X_{11}|T] = E[lambda|T] + Var[lambda|T]  (Law of Total Variance)\n")
cat(sprintf("               = %.4f    + %.4f     = %.4f\n",
            E_lambda_post, Var_lambda_post, E_lambda_post + Var_lambda_post))
cat("  The classical approach ignores parameter uncertainty (Var[lambda|T] term).\n\n")

# --- P(X_{11} > 5 | T) ---
prob_gt5_bayes <- 1 - pnbinom(5, size = r_nb, prob = p_nb)
prob_gt5_class <- 1 - ppois(5, lambda = lambda_mle)

cat("=== P(X_{11} > 5 | T) ===\n")
cat(sprintf("Bayesian (NegBin): %.4f\n", prob_gt5_bayes))
cat(sprintf("Classical (Poisson): %.4f\n", prob_gt5_class))
cat("Classical gives higher probability because its mean (3.5) > Bayesian mean (3.364).\n")
cat("The prior Gamma(2,1) shrinks the Bayesian estimate toward the prior mean (2).\n")

# =============================================================================
# Plot: Prior vs Posterior of lambda
# =============================================================================

lambda_grid <- seq(0, 8, length.out = 500)

plot(lambda_grid, dgamma(lambda_grid, shape = 2, rate = 1),
     type = "l", col = "orange", lwd = 2,
     main = "Prior vs Posterior Distribution of lambda",
     xlab = expression(lambda), ylab = "Density",
     ylim = c(0, 1.6))
lines(lambda_grid, dgamma(lambda_grid, shape = alpha_post, rate = beta_post),
      col = "steelblue", lwd = 2)
abline(v = 2,              col = "orange",   lty = 2, lwd = 1.5)  # prior mean / H0
abline(v = alpha_post / beta_post, col = "steelblue", lty = 2, lwd = 1.5)  # posterior mean
abline(v = T_obs / n_obs,  col = "darkgreen", lty = 2, lwd = 1.5)  # MLE
legend("topright",
       legend = c(sprintf("Prior: Gamma(2,1), mean=2"),
                  sprintf("Posterior: Gamma(%d,%d), mean=%.2f", alpha_post, beta_post,
                          alpha_post / beta_post),
                  sprintf("MLE = %.1f", T_obs / n_obs)),
       col  = c("orange", "steelblue", "darkgreen"),
       lty  = c(1, 1, 2), lwd = c(2, 2, 1.5))
