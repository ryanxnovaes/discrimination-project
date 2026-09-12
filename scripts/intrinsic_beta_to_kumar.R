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


# ============================================================
# Beta expectation involving the Kumaraswamy shape parameter
# ============================================================
#
# Unlike the previous two expectations, this quantity does not
# have a simple closed-form expression for general p.
#
# E_Beta[log(1 - X^p)] = integral_0^1 log(1 - x^p) f_Beta(x) dx.
#
# Direct numerical evaluation of this integral can be unstable
# when the Beta density is highly concentrated near 0 or 1,
# because the density may be singular at the boundaries.
#
# To improve numerical stability, use the decomposition
#
# log(1 - x^p) = log(p) + log(1 - x) + log{(1 - x^p) / [p(1 - x)]}.
#
# Therefore,
#
# E_Beta[log(1 - X^p)] = log(p) + E_Beta[log(1 - X)] + E_Beta[r_p(X)],
#
# where
#
# r_p(x) = log{(1 - x^p) / [p(1 - x)]}.
#
# The expectation E_Beta[log(1 - X)] is available analytically.
# Only the bounded remainder r_p(X) is evaluated numerically.
#
# For the numerical component, the probability integral
# transform is used:
#
# U ~ Uniform(0, 1),  X = Q_Beta(U),
#
# so that
#
# E_Beta[r_p(X)] = integral_0^1 r_p(Q_Beta(u)) du.
#
# This avoids direct integration against a potentially singular
# Beta density. The boundary limits of the remainder are finite:
#
# r_p(0) = -log(p),    r_p(1) = 0.
#
# These limits are explicitly used when the numerical Beta
# quantile is returned as exactly 0 or 1.
# ============================================================

beta_expected_log1m_xp <- function(p, mu, phi, rel.tol = 1e-10) {
  
  if (!is.finite(p) || p <= 0)
    stop("p must be positive")
  
  # Analytic component: E_Beta[log(1 - X)].
  elog_1m_x <- beta_expected_log1m_x(mu = mu, phi = phi)
  
  # Bounded remainder evaluated on the probability scale.
  remainder <- function(u) {
    
    x <- q_beta(p = u, mu = mu, phi = phi)
    out <- numeric(length(x))
    
    at_zero <- x <= 0
    at_one <- x >= 1
    interior <- !at_zero & !at_one
    
    # Use the analytical boundary limits of r_p(x).
    out[at_zero] <- -log(p)
    out[at_one] <- 0
    
    if (any(interior)) {
      
      log_x <- log(x[interior])
      
      out[interior] <- log1mexp(p * log_x) - log1mexp(log_x) - log(p)
    }
    
    out
  }
  
  remainder_expectation <- integrate_unit_interval(f = remainder, 
                                                   rel.tol = rel.tol)
  
  log(p) + elog_1m_x + remainder_expectation
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
# For each p, q is replaced by its analytical profile solution
# q*(p). The KL divergence therefore becomes a function of p
# alone and is optimized on the log scale, p = exp(log_p).
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
# Compute the KL projection of Beta(mu, phi) onto the
# Kumaraswamy family.
#
# Since q*(p) is available analytically, the projection reduces
# to one-dimensional optimization over log(p). The initial
# interval p in [0.01, 100] is only a numerical search region
# and is expanded whenever the optimum lies near a boundary.
#
# The final interval and number of expansions are returned as
# numerical diagnostics.
# ============================================================

project_beta_to_kumar <- function(mu, phi, log_p_interval = NULL,
                                  rel.tol = 1e-10, opt.tol = 1e-12,
                                  boundary_margin = 0.25,
                                  expansion = 2,
                                  max_expansions = 8) {
  
  # Initial search interval on the log(p) scale
  if (is.null(log_p_interval)) {
    log_p_interval <- log(c(0.01, 100))
  }
  
  interval <- log_p_interval
  
  expansions <- 0L
  
  repeat {
    
    opt <- optimize(
      f = beta_to_kumar_kl_profile,
      interval = interval,
      mu = mu,
      phi = phi,
      rel.tol = rel.tol,
      tol = opt.tol
    )
    
    log_p_star <- unname(opt$minimum)
    
    # Check whether optimum is close to a boundary
    near_left <- (log_p_star - interval[1]) < boundary_margin
    near_right <- (interval[2] - log_p_star) < boundary_margin
    
    
    if (!near_left && !near_right)
      break
    
    if (expansions >= max_expansions)
      break
    
    if (near_left)
      interval[1] <- interval[1] - expansion
    
    if (near_right)
      interval[2] <- interval[2] + expansion
    
    expansions <- expansions + 1L
  }
  
  # Recover the optimal Kumaraswamy shape parameters.
  p_star <- exp(log_p_star)
  q_star <- beta_to_kumar_q_star(p = p_star, mu = mu, phi = phi,
                                 rel.tol = rel.tol)
  
  # Convert the projection to the median-dispersion
  # parameterization used for reporting.
  pars_k <- kumar_from_pq(p = p_star, q = q_star)
  
  list(
    p = unname(p_star),
    q = unname(q_star),
    omega = unname(pars_k["omega"]),
    dp = unname(pars_k["dp"]),
    D_star = unname(opt$objective),
    log_p = unname(log_p_star),
    search_interval = unname(interval),
    expansions = expansions,
    boundary_hit = near_left || near_right,
    optim = opt
  )
}


# ============================================================
# Complete intrinsic separation: Beta -> Kumaraswamy
# ============================================================
#
# Compute the separation between a Beta distribution and its
# KL-optimal Kumaraswamy projection.
#
# The function combines the KL divergence, Hellinger distance,
# overlap coefficient, log-likelihood-ratio variance, and the
# standardized separation delta = D_star / sqrt(v).
# ============================================================
intrinsic_beta_to_kumar <- function(mu, phi, log_p_interval = NULL,
                                    rel.tol = 1e-10, zero.tol = 1e-10) {
  
  # Obtain the KL-optimal Kumaraswamy projection.
  proj <- project_beta_to_kumar(
    mu = mu,
    phi = phi,
    log_p_interval = log_p_interval,
    rel.tol = rel.tol
  )
  
  # Define the true Beta density and its log-density.
  d_true <- function(x) {
  
    d_beta(x = x, mu = mu, phi = phi)
  }
  
  log_d_true <- function(x) {
    
    d_beta(x = x, mu = mu, phi = phi, log = TRUE)
  }
  
  # Define the projected Kumaraswamy density and its log-density.
  d_comp <- function(x) {
    
    d_kumar(x = x, omega = proj$omega, dp = proj$dp
    )
  }
  
  log_d_comp <- function(x) {
    
    d_kumar(x = x, omega = proj$omega, dp = proj$dp, log = TRUE)
  }
  
  # Compute the Hellinger distance between the true and projected densities.
  H_star <- hellinger_distance(log_d_true = log_d_true, log_d_comp = log_d_comp,
                               rel.tol = rel.tol
  )
  
  # Compute the overlap coefficient between the true and projected densities.
  OVL_star <- overlap_coefficient(log_d_true = log_d_true, log_d_comp = log_d_comp,
                                  rel.tol = rel.tol
  )
  
  # KL divergence at the optimal Kumaraswamy projection.
  D_star <- unname(proj$D_star)
  
  # Variance of the log-likelihood ratio under the Beta truth.
  v <- loglik_ratio_variance(
    d_true = d_true,
    log_d_true = log_d_true,
    log_d_comp = log_d_comp,
    D_star = D_star,
    rel.tol = rel.tol,
    zero.tol = zero.tol
  )
  
  # Treat sufficiently small numerical values as zero.
  if (is.finite(D_star) && abs(D_star) < zero.tol) {
    D_star <- 0
  }
  
  if (is.finite(v) && abs(v) < zero.tol) {
    v <- 0
  }
  
  # Standardized separation based on the KL divergence and
  # the variance of the log-likelihood ratio.
  delta <- if (is.finite(D_star) && is.finite(v) && D_star > zero.tol && v > zero.tol) {
    D_star / sqrt(v)
  } else {
    NA_real_
  }
  
  # Return the projection parameters and all separation measures.
  list(
    truth = "Beta",
    mu = unname(mu),
    phi = unname(phi),
    p_star = unname(proj$p),
    q_star = unname(proj$q),
    omega_star = unname(proj$omega),
    dp_star = unname(proj$dp),
    D_star = unname(D_star),
    H_star = unname(H_star),
    OVL_star = unname(OVL_star),
    v = unname(v),
    delta = unname(delta)
  )
}