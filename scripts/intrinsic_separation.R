# ============================================================
# Intrinsic separation between Beta and Kumaraswamy models
# ============================================================

source("scripts/distributions.R")
source("scripts/estimation.R")

# ============================================================
#  Numerical integration on (0, 1)
# ============================================================

# Splits (0,1) into two subintervals to improve numerical 
# stability near the boundaries x = 0 and x = 1.

integrate_unit_interval <- function(f, rel.tol = 1e-10, subdivisions = 2000L) {
  
  # Maximum number of adaptive subdivisions allowed by integrate().
  left <- integrate(
    f = f,
    lower = 0,
    upper = 0.5,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )
  
  # Maximum number of adaptive subdivisions allowed by integrate().
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