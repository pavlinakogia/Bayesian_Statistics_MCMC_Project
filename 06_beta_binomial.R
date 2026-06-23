# =============================================================================
# 06_beta_binomial.R
# Bayesian Statistics and MCMC - Exercise 3
#
# Bayesian Beta-Binomial analysis of esophageal cancer risk by alcohol
# consumption category (Breslow & Day, 1980 -- built-in R dataset 'esoph').
#
# Focus: age group 45-54 years.
# Response:     ncases  ~ Binomial(n_i, p_i)
# Explanatory:  alcgp (4 levels: 0-39, 40-79, 80-119, 120+ g/day)
#
# Priors:
#   p[1] ~ Beta(2, 38)   -- informative: based on hypothetical past study
#                           with 2 cases in 40 trials => prior mean = 0.05
#                           (low cancer probability for low alcohol consumption)
#   p[2], p[3], p[4] ~ Beta(1, 1) = Uniform(0,1)  -- diffuse (no prior info)
#
# Conjugate update: Beta(a,b) + Binomial(n,p) => Beta(a+y, b+n-y)
#
# Odds Ratios (vs reference category 0-39 g/day) are computed from MCMC samples.
# This script provides the analytical posterior update and simulation-based OR.
# A WinBUGS model is provided in models/winbugs_exercise3.txt for full MCMC.
# =============================================================================

# =============================================================================
# Data: esoph dataset, age group 45-54
# =============================================================================

data(esoph)
esoph_45_54 <- subset(esoph, agegp == "45-54")

# Aggregate by alcohol consumption category
library(dplyr)
agg <- esoph_45_54 %>%
  group_by(alcgp) %>%
  summarise(cases  = sum(ncases),
            total  = sum(ncases + ncontrols),
            .groups = "drop")

cat("=== Data: Esophageal Cancer by Alcohol Consumption (Age 45-54) ===\n")
print(agg)
cat("\n")

y <- agg$cases    # number of cancer cases per category
n <- agg$total    # total individuals per category
N <- nrow(agg)    # number of categories = 4

# Observed proportions (= MLE)
p_mle <- y / n
cat("Observed proportions (MLE):\n")
for (i in 1:N) cat(sprintf("  %s: %d/%d = %.4f\n", agg$alcgp[i], y[i], n[i], p_mle[i]))

# =============================================================================
# Bayesian Update (analytical -- Beta is conjugate to Binomial)
# =============================================================================

# Prior parameters
alpha_prior <- c(2, 1, 1, 1)   # Beta(2,38) for dose=1; Beta(1,1) for others
beta_prior  <- c(38, 1, 1, 1)

# Posterior parameters: Beta(alpha + y, beta + n - y)
alpha_post <- alpha_prior + y
beta_post  <- beta_prior  + (n - y)

# Posterior means and 95% CrI
p_bayes_mean <- alpha_post / (alpha_post + beta_post)
p_bayes_ci   <- matrix(NA, nrow = N, ncol = 2)
for (i in 1:N) {
  p_bayes_ci[i, ] <- qbeta(c(0.025, 0.975), alpha_post[i], beta_post[i])
}

cat("\n=== Posterior Parameters ===\n")
for (i in 1:N) {
  cat(sprintf("  %s: Beta(%d, %d)  =>  mean = %.4f, 95%% CrI = (%.4f, %.4f)\n",
              agg$alcgp[i], alpha_post[i], beta_post[i],
              p_bayes_mean[i], p_bayes_ci[i, 1], p_bayes_ci[i, 2]))
}

# =============================================================================
# Comparison Table: MLE vs Bayesian
# =============================================================================

results <- data.frame(
  Category    = agg$alcgp,
  Cases       = y,
  Total       = n,
  MLE         = round(p_mle, 4),
  Bayes_Mean  = round(p_bayes_mean, 4),
  CrI_lower   = round(p_bayes_ci[, 1], 4),
  CrI_upper   = round(p_bayes_ci[, 2], 4)
)

cat("\n=== Comparison: MLE vs Bayesian ===\n")
print(results, row.names = FALSE)

cat("\nKey observation: For dose=1 (0-39 g/day), Bayesian mean (", round(p_bayes_mean[1], 4),
    ") > MLE (", round(p_mle[1], 4), ")\n")
cat("This is due to the informative prior Beta(2,38) with mean 0.05,\n")
cat("which pulls the estimate upward from the MLE (1/78 = 0.013).\n")
cat("With only 1 case in 78 trials, prior information has large influence.\n")

# =============================================================================
# Odds Ratios via Monte Carlo simulation
# (Same result as WinBUGS -- analytical Beta posterior sampled directly)
# =============================================================================

set.seed(42)
n_sim <- 50000

# Sample from posterior distributions
p_sim <- matrix(NA, nrow = n_sim, ncol = N)
for (i in 1:N) {
  p_sim[, i] <- rbeta(n_sim, alpha_post[i], beta_post[i])
}

# Compute Odds Ratios vs reference category (dose=1)
odds     <- p_sim / (1 - p_sim)
OR2_sim  <- odds[, 2] / odds[, 1]
OR3_sim  <- odds[, 3] / odds[, 1]
OR4_sim  <- odds[, 4] / odds[, 1]

or_results <- data.frame(
  Comparison = c("40-79 vs 0-39", "80-119 vs 0-39", "120+ vs 0-39"),
  OR_mean    = round(c(mean(OR2_sim), mean(OR3_sim), mean(OR4_sim)), 2),
  CI_2.5     = round(c(quantile(OR2_sim, 0.025), quantile(OR3_sim, 0.025),
                        quantile(OR4_sim, 0.025)), 2),
  CI_97.5    = round(c(quantile(OR2_sim, 0.975), quantile(OR3_sim, 0.975),
                        quantile(OR4_sim, 0.975)), 2)
)

cat("\n=== Bayesian Odds Ratios (vs 0-39 g/day reference) ===\n")
print(or_results, row.names = FALSE)
cat("\nAll 95% CrIs exclude 1 => statistically significant association\n")
cat("between alcohol consumption and esophageal cancer for all categories.\n")
cat("\nNote: Wide CrI for 120+ vs 0-39 reflects:\n")
cat("  (1) Small sample size in 120+ category (n=15)\n")
cat("  (2) p4 near 1 and p1 near 0 => ratio p/(1-p) is non-linearly sensitive\n")
cat("  (3) Two uncertain quantities in the ratio => uncertainty multiplies\n")

# =============================================================================
# Plots
# =============================================================================

# Plot 1: Prior vs Posterior for p1 (illustrates informative prior effect)
lambda_grid <- seq(0, 0.2, length.out = 500)
plot(lambda_grid, dbeta(lambda_grid, alpha_prior[1], beta_prior[1]),
     type = "l", col = "orange", lwd = 2,
     main = "Prior vs Posterior: p1 (0-39 g/day)",
     xlab = "p1", ylab = "Density")
lines(lambda_grid, dbeta(lambda_grid, alpha_post[1], beta_post[1]),
      col = "steelblue", lwd = 2)
abline(v = p_mle[1],        col = "darkgreen", lty = 2, lwd = 1.5)
abline(v = p_bayes_mean[1], col = "steelblue", lty = 2, lwd = 1.5)
legend("topright",
       legend = c("Prior Beta(2,38)", "Posterior Beta(3,115)",
                  sprintf("MLE = %.3f", p_mle[1]),
                  sprintf("Bayes mean = %.3f", p_bayes_mean[1])),
       col  = c("orange", "steelblue", "darkgreen", "steelblue"),
       lty  = c(1, 1, 2, 2), lwd = c(2, 2, 1.5, 1.5))

# Plot 2: Comparison Bayesian vs MLE across categories
plot(1:N, p_bayes_mean, pch = 17, col = "steelblue", cex = 1.5,
     ylim = c(0, 1.05),
     xaxt = "n",
     main = "Probability of Esophageal Cancer by Alcohol Consumption\n(Age 45-54)",
     xlab = "Alcohol Consumption Category",
     ylab = "Probability of Cancer")
axis(1, at = 1:N, labels = agg$alcgp)
segments(1:N, p_bayes_ci[, 1], 1:N, p_bayes_ci[, 2],
         col = "steelblue", lwd = 2)
points(1:N, p_mle, pch = 4, col = "red", cex = 1.5)
points(1:N, y / n, pch = 3, col = "darkgreen", cex = 1.2)
legend("topleft",
       legend = c("Bayesian mean + 95% CrI", "MLE", "Observed proportion"),
       pch    = c(17, 4, 3),
       col    = c("steelblue", "red", "darkgreen"))
