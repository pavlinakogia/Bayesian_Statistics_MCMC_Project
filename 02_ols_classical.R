# =============================================================================
# 02_ols_classical.R
# Bayesian Statistics and MCMC - Exercise 1
#
# Classical OLS estimation for comparison with Bayesian results.
# The OLS estimator beta_hat = (X'X)^{-1} X'y is also the MLE under
# Normal errors -- this equivalence holds specifically because the
# likelihood's exponent contains ||y - X*beta||^2.
# =============================================================================

# Run data simulation first
source("01_data_simulation.R")

# --- Design matrix with intercept column ---
X_full <- cbind(1, X)
colnames(X_full) <- c("intercept", paste0("X", 1:15))

# --- OLS estimation via lm() ---
# lm() internally uses QR decomposition (numerically stable)
# equivalent to beta_hat = (X'X)^{-1} X'y
df <- as.data.frame(X_full)
df$y <- y
fit_ols <- lm(y ~ . - 1, data = df)  # -1 because intercept already in X_full

# --- Extract results ---
ols_coef    <- coef(fit_ols)
ols_se      <- summary(fit_ols)$coefficients[, "Std. Error"]
ols_ci      <- confint(fit_ols, level = 0.95)
ols_sigma2  <- sum(residuals(fit_ols)^2) / fit_ols$df.residual  # s^2 / (n-p)
ols_s2_raw  <- sum(residuals(fit_ols)^2)  # raw RSS = s^2 in the notation of the assignment

# --- Summary table ---
results_ols <- data.frame(
  True      = true_beta,
  OLS       = round(ols_coef, 4),
  CI_lower  = round(ols_ci[, 1], 4),
  CI_upper  = round(ols_ci[, 2], 4),
  row.names = names(ols_coef)
)

cat("=== OLS Results ===\n")
print(results_ols)
cat("\nOLS estimate of sigma^2:", round(ols_sigma2, 4), "\n")
cat("True sigma^2:", true_sigma2, "\n")
cat("RSS (s^2):", round(ols_s2_raw, 4), "\n")

# --- Variable selection: keep variables whose 95% CI excludes zero ---
keep_ols <- rownames(ols_ci)[ols_ci[, 1] > 0 | ols_ci[, 2] < 0]
cat("\nVariables with CI excluding zero (OLS):", paste(keep_ols, collapse = ", "), "\n")
