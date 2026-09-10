source("scripts/distributions.R")

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
# Hellinger distance
# ============================================================

# Compute the Hellinger distance between two densities on (0, 1).
#
# H(f, g) = sqrt(1 - integral_0^1 sqrt(f(x) g(x)) dx).
#
# The densities are supplied on the log scale. The integrand is
# therefore evaluated as exp((log f(x) + log g(x)) / 2).
hellinger_distance <- function(log_d_true, log_d_comp, rel.tol = 1e-10) {
  
  integrand <- function(x) {
    
    log_f <- log_d_true(x)
    log_g <- log_d_comp(x)
    
    exp(0.5 * (log_f + log_g))
  }
  
  # Compute the Hellinger affinity by numerical integration.
  affinity <- integrate_unit_interval(f = integrand, rel.tol = rel.tol)
  
  # Correct small numerical violations of the theoretical upper bound.
  if (affinity > 1 && affinity - 1 < 1e-8)
    affinity <- 1
  
  if (affinity < 0 || affinity > 1 + 1e-8)
    warning("Invalid Hellinger affinity")
  
  # Hellinger distance from the affinity.
  h2 <- max(0, 1 - affinity)
  
  sqrt(h2)
}

# - hellinger_distance
# - overlap_coefficient
# - loglik_ratio_variance