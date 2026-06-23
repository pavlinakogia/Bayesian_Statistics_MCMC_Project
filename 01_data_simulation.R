# =============================================================================
# 01_data_simulation.R
# Bayesian Statistics and MCMC - Exercise 1
# 
# Simulates data for a multiple linear regression problem:
#   Y = X*beta + epsilon, epsilon ~ N(0, sigma^2)
#
# Design:
#   - n = 50 observations, p = 15 predictors + intercept
#   - X1-X10: independent standard Normal
#   - X11-X15: linear functions of X1-X5 + noise (induces multicollinearity)
#   - True coefficients: beta0=4, beta1=2, beta5=-1, beta7=1.5, beta11=1, beta13=0.5
#   - True sigma^2 = 6.25 (sigma = 2.5)
# =============================================================================

set.seed(123)
n <- 50

# --- Simulate X1-X10: independent N(0,1) ---
X1 <- matrix(rnorm(n * 10), nrow = n, ncol = 10)

# --- Simulate X11-X15: linear combinations of X1-X5 + noise ---
# This induces multicollinearity between X1-X5 and X11-X15
X2 <- matrix(NA, nrow = n, ncol = 5)
for (j in 1:5) {
  X2[, j] <- rnorm(n,
                   mean = 0.3 * X1[, 1] + 0.5 * X1[, 2] +
                     0.7 * X1[, 3] + 0.9 * X1[, 4] + 1.5 * X1[, 5],
                   sd = 1)
}

# --- Full design matrix (without intercept column) ---
X <- cbind(X1, X2)

# --- Simulate response variable ---
# True model: Y = 4 + 2*X1 - X5 + 1.5*X7 + X11 + 0.5*X13 + epsilon
# Note: X11 = X2[,1], X13 = X2[,3] in the full X matrix
y <- rnorm(n,
           mean = 4 + 2 * X1[, 1] - X1[, 5] + 1.5 * X1[, 7] +
             X2[, 1] + 0.5 * X2[, 3],
           sd = 2.5)

# --- True coefficient vector (for reference) ---
true_beta <- c(4,    # beta0 (intercept)
               2,    # beta1
               0,    # beta2
               0,    # beta3
               0,    # beta4
               -1,   # beta5
               0,    # beta6
               1.5,  # beta7
               0,    # beta8
               0,    # beta9
               0,    # beta10
               1,    # beta11
               0,    # beta12
               0.5,  # beta13
               0,    # beta14
               0)    # beta15

true_sigma2 <- 6.25

cat("Data simulation complete.\n")
cat("n =", n, "observations,", ncol(X), "predictors\n")
cat("True sigma^2 =", true_sigma2, "\n")
cat("Non-zero coefficients: beta0, beta1, beta5, beta7, beta11, beta13\n")
