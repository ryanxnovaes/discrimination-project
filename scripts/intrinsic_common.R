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

# - hellinger_distance
# - overlap_coefficient
# - loglik_ratio_variance