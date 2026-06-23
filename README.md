# Bayesian Statistics and MCMC

**MSc in Mathematical Modelling in Modern Technologies and Finance**  
National Technical University of Athens (NTUA)  
Spring Semester 2026

---

## Overview

This project implements Bayesian statistical methods and MCMC algorithms from scratch in R, covering four interconnected exercises:

1. **Bayesian Linear Regression** — conjugate Normal-Inverse Gamma model, Gibbs sampler
2. **Poisson-Gamma Model** — Bayesian hypothesis testing with Bayes Factors
3. **Beta-Binomial Model** — informative priors, Odds Ratios (esophageal cancer data)
4. **Metropolis-Hastings** — custom MH algorithm under a Jeffreys prior

The emphasis is on understanding the theory (deriving posteriors analytically) and then verifying results computationally through MCMC.

---

## Key Skills Demonstrated

- **From-scratch MCMC**: Gibbs sampler and Metropolis-Hastings implemented in pure R without MCMC packages
- **Bayesian model comparison**: Bayes Factors, posterior odds, Jeffreys scale interpretation
- **Conjugate analysis**: Normal-IG, Gamma-Poisson, Beta-Binomial closed-form posteriors
- **MCMC diagnostics**: trace plots, ACF, ergodic mean plots, running quantiles
- **Comparison with classical statistics**: OLS vs Bayesian estimates, shrinkage effect

---

## Repository Structure

```
├── 01_data_simulation.R       # Simulate regression data (n=50, p=15, multicollinearity)
├── 02_ols_classical.R         # Classical OLS estimation (baseline comparison)
├── 03_gibbs_sampler.R         # Custom Gibbs sampler for Bayesian linear regression ★
├── 04_metropolis_hastings.R   # Custom MH algorithm under Jeffreys prior ★
├── 05_bayes_factor.R          # Bayesian hypothesis test + posterior predictive
├── 06_beta_binomial.R         # Beta-Binomial analysis (esophageal cancer dataset)
├── models/
│   ├── winbugs_exercise1.txt  # WinBUGS model: Bayesian linear regression
│   └── winbugs_exercise3.txt  # WinBUGS model: Beta-Binomial
└── Kogia-Pavlina.pdf          # Full assignment report (with all derivations in LaTeX)
```

★ = most technically interesting files (MCMC from scratch)

---

## Exercises in Detail

### Exercise 1 — Bayesian Linear Regression

**Model**: `Y = Xβ + ε`, `ε ~ N(0, σ²Iₙ)`, `n=50`, `p=15`

**Priors** (conjugate Normal-Inverse Gamma):
- `β | σ² ~ N₁₆(0, σ²M⁻¹)` with `M⁻¹ = 10000·I₁₆` (diffuse)
- `σ² ~ IG(0.0001, 0.0001)`

**Key results**:
- Derived analytically that `OLS = MLE` under Normal errors
- Posterior `β | σ², y, X ~ N₁₆(β*, σ²(M + X'X)⁻¹)` where `β*` is a weighted average of OLS and prior mean
- Posterior `σ² | y, X ~ IG(n/2 + a, b + s²/2 + C/2)`
- Bayesian posterior mean ≈ OLS (diffuse prior → data dominate)
- **Bayesian σ² estimate < OLS**: denominator is `n ≈ 48` vs `n-p = 34` (marginalization effect)
- Simulated data includes intentional multicollinearity (X11–X15 are linear functions of X1–X5), causing high posterior variance for related coefficients

**`03_gibbs_sampler.R`** — Two-block Gibbs sampler alternating between:
1. `β | σ², y, X` → `rmvnorm()` (multivariate Normal)
2. `σ² | β, y, X` → `1/rgamma()` (Inverse Gamma via reciprocal)

---

### Exercise 2 — Bayesian Hypothesis Testing (Poisson)

**Data**: `T = 35` errors observed over `n = 10` days, `T ~ Poisson(nλ)`

**Test**: `H₀: λ = 2` vs `H₁: λ ≠ 2`

**Prior under H₁**: `λ ~ Gamma(2, 1)` — weak belief that λ ≈ 2 (same mean as H₀, large variance)

**Results**:
- `BF₁₀ = 15.45` → strong evidence for H₁ (Jeffreys scale: 10–30)
- `PO₁₀ = BF₁₀ = 15.45` (equal prior probabilities → posterior odds = Bayes Factor)
- Posterior predictive: `X₁₁ | T ~ NegBin(37, 11/12)` (Poisson-Gamma mixture)
- `Var[X₁₁|T] = E[λ|T] + Var[λ|T] = 3.364 + 0.306 = 3.669 > 3.5` (classical)
- **Classical variance underestimates** by ignoring parameter uncertainty (Law of Total Variance)

---

### Exercise 3 — Beta-Binomial (Esophageal Cancer)

**Data**: `esoph` dataset (R built-in), age group 45–54, Breslow & Day (1980)

**Model**: `Yᵢ ~ Binomial(nᵢ, pᵢ)` for 4 alcohol consumption categories

**Priors**:
- `p₁ ~ Beta(2, 38)`: informative (hypothetical past study — 2 cases in 40 trials, prior mean = 0.05)
- `p₂, p₃, p₄ ~ Beta(1, 1)`: diffuse

**Key results**:

| Category | MLE | Bayes Mean | OR (vs 0–39) |
|---|---|---|---|
| 0–39 g/day | 0.013 | 0.025 | — |
| 40–79 g/day | 0.247 | 0.253 | 19.88 |
| 80–119 g/day | 0.308 | 0.317 | 27.60 |
| 120+ g/day | 0.867 | 0.824 | 403.78 |

- Largest Bayesian–classical difference at dose=1: informative prior pulls estimate from 0.013 → 0.025 (only 1 case in 78 trials — prior dominates)
- Wide CrI for OR(120+): small sample (n=15) + p₄ near 1 + p₁ near 0 → uncertainty multiplies in ratio

---

### Exercise 4 — Metropolis-Hastings under Jeffreys Prior

**Jeffreys prior** for `λ`: `p(λ) ∝ λ^{-1/2}` (improper, but posterior is proper)

Derived via Fisher Information: `I(λ) = E[-∂²/∂λ² log f(x|λ)] = 1/λ`

**Analytical posterior**: `λ | T ~ Gamma(35.5, 10)`

**MH algorithm** (`04_metropolis_hastings.R`):
- Proposal: Normal random walk `λ* ~ N(λ_current, σ² = 8)`
- Acceptance rate ≈ 25% (as specified)
- Symmetric proposal → acceptance ratio simplifies to posterior ratio (no proposal terms)
- All computation on log scale for numerical stability

**Results**: MH samples match `Gamma(35.5, 10)` closely (mean: 3.527 vs 3.550 analytical)

**95% Credible Interval**: `(2.502, 4.729)` — H₀ value `λ=2` falls outside → consistent with Exercise 2

---

## How to Run

### Requirements

```r
install.packages(c("mvtnorm", "dplyr"))
```

> **Note on WinBUGS**: The `.txt` model files in `models/` are provided for reference and reproducibility of the full analysis. WinBUGS runs on Windows only. The core MCMC algorithms (Gibbs sampler and Metropolis-Hastings) are fully implemented in pure R in `03_gibbs_sampler.R` and `04_metropolis_hastings.R` — no WinBUGS required to run these.

### Run order

```r
source("01_data_simulation.R")    # must run first (generates X and y)
source("02_ols_classical.R")      # depends on 01
source("03_gibbs_sampler.R")      # depends on 01 and 02
source("04_metropolis_hastings.R") # standalone
source("05_bayes_factor.R")        # standalone
source("06_beta_binomial.R")       # standalone
```

---

## Theoretical Highlights

### Why OLS = MLE under Normal errors

Maximizing `ℓ(β, σ²|y, X) ∝ exp(-||y-Xβ||²/2σ²)` over β is equivalent to minimizing `||y-Xβ||²` — the OLS criterion. This equivalence is specific to Normal errors; under Laplace errors, MLE minimizes sum of absolute residuals instead.

### Completing the square

The key algebraic step in deriving the posteriors:

```
||y - Xβ||² = s² + (β - β̂)'(X'X)(β - β̂)
```

where `s² = ||y - Xβ̂||²` (RSS) and the cross term vanishes by the normal equations `X'(y - Xβ̂) = 0` (the OLS residual is orthogonal to the column space of X).

### Why Bayesian σ² < OLS σ²

OLS: `s²/(n-p) = RSS/34`  
Bayesian posterior mean: `≈ RSS/48`  
Ratio: `34/48 ≈ 0.708` — confirmed empirically (3.87/5.47 ≈ 0.708)

The proper (even if diffuse) prior on β changes the effective degrees of freedom from `n-p` to approximately `n` when marginalizing over β.

### Gibbs vs Metropolis-Hastings

| | Gibbs | Metropolis-Hastings |
|---|---|---|
| Requires | Analytical full conditionals | Only unnormalized posterior |
| Accept/reject | Never (always accepted) | Yes (25% accepted here) |
| Autocorrelation | Near-zero from lag 1 | Significant (~10 lags) |
| Effective sample size | ≈ n_iter | ≪ n_iter |
| Applicable when | Conjugate priors known | Any posterior |

---

## Report

The full derivations (likelihood factorization, completing the square, posterior updates, Bayes Factor computation, Jeffreys prior derivation) are in `Kogia-Pavlina.pdf`, typeset in LaTeX.

---

*MSc Mathematical Modelling — NTUA, June 2026*
