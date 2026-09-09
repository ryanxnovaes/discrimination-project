# ============================================================
# Intrinsic separation between Beta and Kumaraswamy models
# ============================================================

source("scripts/distributions.R")
source("scripts/intrinsic_common.R")

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


# ============================================================
# Beta -> Kumaraswamy:
# Profiled Kumaraswamy parameter q
# ============================================================
#
# Model: X ~ Beta(alpha, beta)
# approximated by a Kumaraswamy distribution with density:
# f_K(x | p, q) = p q x^(p - 1) (1 - x^p)^(q - 1).
#
# For fixed p, the expected Kumaraswamy log-density under
# Beta truth depends on q through
# log(q) + (q - 1) E_Beta[log(1 - X^p)].
#
# Differentiating the expected log-density with respect to q
# and setting the derivative equal to zero gives the
# closed-form profile solution
#
#   q*(p) = -1 / E_Beta[log(1 - X^p)].
#
# Therefore, q is obtained analytically for each fixed p.
# This is a profile optimization: rather than optimizing the
# KL divergence jointly over the two parameters (p, q), the
# optimal q is first expressed as a function of p.
#
# The original two-dimensional optimization is thus reduced
# to a one-dimensional optimization over p:
#
#   (p, q)  -->  q*(p)
#           -->  D_KL(Beta || Kumaraswamy)
#           -->  optimize only over p.
#
# This avoids a numerical optimization over q and exploits
# the analytical solution available for fixed p.
# ============================================================

beta_to_kumar_q_star <- function(p, mu, phi, rel.tol = 1e-10) {
  
  elog_1m_xp <- beta_expected_log1m_xp(p = p, mu = mu, phi = phi, 
                                       rel.tol = rel.tol)
  
  q_star <- -1 / elog_1m_xp
  
  if (!is.finite(q_star) || q_star <= 0)
    stop("Invalid profiled Kumaraswamy q")
  
  q_star
}


# Expected Kumaraswamy log-density under Beta truth
beta_expected_log_kumar <- function(p, mu, phi, rel.tol = 1e-10) {
  
  if (!is.finite(p) || p <= 0)
    return(-Inf)
  
  elog_x <- beta_expected_log_x(mu = mu, phi = phi)
  elog_1m_xp <- beta_expected_log1m_xp(p = p, mu = mu, phi = phi, 
                                       rel.tol = rel.tol)
  
  # For fixed p, obtain the KL-optimal profiled value q*(p).
  q_star <- -1 / elog_1m_xp
  
  if (!is.finite(q_star) || q_star <= 0)
    return(-Inf)
  
  # Evaluate the expected Kumaraswamy log-density under Beta truth.
  log(p) + log(q_star) + (p - 1) * elog_x + (q_star - 1) * elog_1m_xp
}


# ============================================================
# Profiled KL divergence: Beta -> Kumaraswamy
# ============================================================
#
# The Kullback-Leibler divergence from the Beta truth to the
# Kumaraswamy approximation is
# D_KL(Beta || Kumaraswamy) = E_Beta[log f_Beta(X)] - E_Beta[log f_K(X)].
#
# For each fixed p, the Kumaraswamy parameter q is replaced by
# its KL-optimal profiled value q*(p). Therefore, the objective
# function depends only on p:
# D_KL^profiled(p) = E_Beta[log f_Beta(X)] - E_Beta[log f_K(X | p, q*(p))].
#
# The profiled KL objective is evaluated on the log scale:
# p = exp(log_p).
# ============================================================

beta_to_kumar_kl_profile <- function(log_p, mu, phi, rel.tol = 1e-10) {
  
  p <- exp(log_p)
  
  elog_beta <- beta_expected_log_beta(mu = mu, phi = phi)
  
  elog_kumar <- beta_expected_log_kumar(p = p, mu = mu, phi = phi, 
                                        rel.tol = rel.tol)
  
  kl <- elog_beta - elog_kumar
  
  if (!is.finite(kl))
    return(.Machine$double.xmax)
  
  unname(kl)
}

# ============================================================
# KL projection: Beta -> Kumaraswamy
# ============================================================
#
# Numerically minimizes the profiled KL objective over log(p)
# to obtain the KL-optimal Kumaraswamy parameter p*.
#
# The corresponding q* is then recovered analytically from
# the profile solution q*(p*), and the optimal parameters are
# converted from (p*, q*) to the (omega, dp) parameterization.
#
# Returns the KL-optimal parameters and optimization results.
# ============================================================

project_beta_to_kumar <- function(mu, phi, log_p_interval = c(-8, 8), rel.tol = 1e-10) {
  
  opt <- optimize(
    f = beta_to_kumar_kl_profile,
    interval = log_p_interval,
    mu = mu,
    phi = phi,
    rel.tol = rel.tol
  )
  
  p_star <- unname(exp(opt$minimum))
  
  q_star <- beta_to_kumar_q_star(p = p_star, mu = mu, phi = phi,rel.tol = rel.tol)
  
  pars_k <- kumar_from_pq(p = p_star, q = q_star)
  
  list(
    p = p_star,
    q = q_star,
    omega = unname(pars_k["omega"]),
    dp = unname(pars_k["dp"]),
    D_star = unname(opt$objective),
    log_p = unname(opt$minimum),
    optim = opt
  )
}