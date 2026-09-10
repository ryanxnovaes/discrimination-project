# ============================================================
# Intrinsic separation between Kumaraswamy and Beta models
# ============================================================

source("scripts/distributions.R")
source("scripts/intrinsic_common.R")

# ============================================================
# Kumaraswamy truth: expected log terms
# ============================================================
#
# Model: X ~ Kumaraswamy(p, q)
#
# Kumaraswamy density:
# f_K(x | p, q) = p q x^(p - 1) (1 - x^p)^(q - 1).
#
# Log-density:
# log f_K(x | p, q) = log(p) + log(q) + (p - 1) log(x) + (q - 1) log(1 - x^p).
#
# The functions below compute expected logarithmic terms
# under Kumaraswamy truth. These quantities are later used
# to evaluate the expected log-density of the Kumaraswamy
# and Beta models.
#
# For E_K[log(X)], a closed-form expression is available:
#
# E_K[log(X)] = [psi(1) - psi(q + 1)] / p,
# where psi(.) is the digamma function.
#
# This result follows from the transformation Y = X^p, for which Y ~ Beta(1, q).
kumar_expected_log_x <- function(omega, dp) {
  
  pars <- kumar_to_pq(omega = omega, dp = dp)
  
  p <- pars["p"]
  q <- pars["q"]

  (digamma(1) - digamma(q + 1)) / p
}


# E_K[log(1 - X)] = integral_0^1 log(1 - x) f_K(x) dx.
kumar_expected_log1m_x <- function(omega, dp, rel.tol = 1e-10) {
  
  integrand <- function(x) {
    
    density <- d_kumar(x = x, omega = omega, dp = dp)
    
    density * log1p(-x)
  }
  
  integrate_unit_interval(f = integrand, rel.tol = rel.tol)
}


# ============================================================
# Kumaraswamy moments
# ============================================================

# E[X^r] = q * B(1 + r / p, q).
# Computed on the log scale for numerical stability.
kumar_moment <- function(r, omega, dp) {
  
  if (!is.finite(r) || r <= 0)
    stop("r must be positive")
  
  pars <- kumar_to_pq(omega = omega, dp = dp)
  
  p <- pars["p"]
  q <- pars["q"]
  
  exp(log(q) + lbeta(1 + r / p,q))
}

# Mean and variance from the first two moments:
# E[X] = m_1,  E[X^2] = m_2,
# Var(X) = E[X^2] - E[X]^2.
kumar_mean_var <- function(omega, dp) {
  
  mean_x <- kumar_moment(r = 1, omega = omega, dp = dp)
  second_moment <- kumar_moment(r = 2, omega = omega,dp = dp)
  
  var_x <- second_moment - mean_x^2
  
  c(mean = mean_x, variance = var_x)
}


