# ============================================================
# Intrinsic separation between Beta and Kumaraswamy models
# ============================================================

source("scripts/distributions.R")
source("scripts/estimation.R")

# ============================================================
# Numerical integration on (0, 1)
# ============================================================

# Splits (0, 1) into two subintervals to improve numerical
# stability near the boundaries. The subdivisions argument
# controls the maximum number of adaptive subdivisions used
# by integrate().

integrate_unit_interval <- function(f, rel.tol = 1e-10, subdivisions = 2000L) {
  
  # Integrate over the lower half of the unit interval.
  left <- integrate(
    f = f,
    lower = 0,
    upper = 0.5,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )
  
  # Integrate over the upper half of the unit interval.
  right <- integrate(
    f = f,
    lower = 0.5,
    upper = 1,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )
  
  # Combine both halves to obtain the integral over (0,1).
  left$value + right$value
}


# ============================================================
# Beta truth: expected log terms
# ============================================================

# Model: X ~ Beta(alpha, beta)
# 
# Beta density:
# f_Beta(x | alpha, beta) = x^(alpha - 1) * (1 - x)^(beta - 1) / B(alpha, beta).
# 
# Log-density:
# log f_Beta(x | alpha, beta) = -log B(alpha, beta) + (alpha - 1) log(x)
#                               + (beta - 1) log(1 - x).
# 
# 
# Mean-precision parameterization:
# alpha = mu * phi
# beta  = (1 - mu) * phi
#
# Under this parameterization:
# log f_Beta(x | mu, phi) = -log B(mu * phi, (1 - mu) * phi) + (mu * phi - 1) log(x)
#                           + ((1 - mu) * phi - 1) log(1 - x).
#
# The functions below compute expected logarithmic terms under
# Beta truth. These quantities are later used to evaluate the
# expected log-density of the Beta and Kumaraswamy models.

# E_Beta[log(X)] = psi(alpha) - psi(alpha + beta)
# where psi(.) is the digamma function.
beta_expected_log_x <- function(mu, phi) {
  
  pars <- beta_to_ab(mu = mu, phi = phi)
  
  alpha <- pars["alpha"]
  beta  <- pars["beta"]
  
  digamma(alpha) - digamma(alpha + beta)
}


# E_Beta[log(1 - X)] = psi(beta) - psi(alpha + beta)
beta_expected_log1m_x <- function(mu, phi) {
  
  pars <- beta_to_ab(mu = mu, phi = phi)
  
  alpha <- pars["alpha"]
  beta  <- pars["beta"]
  
  digamma(beta) - digamma(alpha + beta)
}


# Unlike the previous two expectations, this quantity is
# evaluated numerically over the unit interval. The parameter
# p is the Kumaraswamy shape parameter appearing in log(1 - X^p).

# E_Beta[log(1 - X^p)] = integral_0^1 log(1 - x^p) f_Beta(x) dx.
beta_expected_log1m_xp <- function(p, mu, phi, rel.tol = 1e-10) {
  
  if (!is.finite(p) || p <= 0)
    stop("p must be positive")
  
  integrand <- function(x) {
    
    log_term <- log1mexp(p * log(x))
    density <- d_beta(x = x, mu = mu, phi = phi)
    
    density * log_term
  }
  
  integrate_unit_interval(f = integrand, rel.tol = rel.tol)
}

# Expected Beta log-density under Beta truth: E_Beta[log f_Beta(X)]
# E_Beta[log f_Beta(X)] = -log B(alpha, beta) + (alpha - 1) E_Beta[log(X)]
#                         + (beta - 1) E_Beta[log(1 - X)]

beta_expected_log_beta <- function(mu, phi) {
  
  pars <- beta_to_ab(mu = mu, phi = phi)
  
  alpha <- pars["alpha"]
  beta  <- pars["beta"]
  
  # Compute E_Beta[log(X)]
  elog_x <- beta_expected_log_x(mu = mu, phi = phi)
  
  # Compute E_Beta[log(1 - X)]
  elog_1mx <- beta_expected_log1m_x(mu = mu, phi = phi)
  
  # Combine the three components of the expected Beta log-density.
  -lbeta(alpha, beta) + (alpha - 1) * elog_x + (beta - 1) * elog_1mx
}