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


# ============================================================
# Starting values for Kumaraswamy -> Beta projection
# ============================================================

# Obtain the Kumaraswamy mean and variance.
kumar_to_beta_start <- function(omega, dp) {
  
  moments <- kumar_mean_var(omega = omega, dp = dp)
  
  mu_start <- unname(moments["mean"])
  var_start <- unname(moments["variance"])
  
  # Recover the Beta precision from its mean-variance relation:
  # phi = mu * (1 - mu) / Var(X) - 1.
  phi_start <- mu_start * (1 - mu_start) / var_start - 1
  
  # Ensure a valid positive precision starting value.
  if (!is.finite(phi_start) || phi_start <= 0) {
    phi_start <- 1
  }
  
  # Convert the mean-precision parameters to Beta shape parameters.
  alpha_start <- mu_start * phi_start
  beta_start  <- (1 - mu_start) * phi_start
  
  # If the moment-based starting values are invalid, use
  # alpha = beta = 1, corresponding to Uniform(0, 1).
  if (!is.finite(alpha_start) || !is.finite(beta_start) || 
      alpha_start <= 0 || beta_start <= 0) {
    alpha_start <- 1
    beta_start <- 1
  }
  
  # Return the shape parameters on the log scale.
  c(log_alpha = log(alpha_start), log_beta = log(beta_start))
}


# ============================================================
# Score equations for Kumaraswamy -> Beta projection
# ============================================================

# The KL projection is obtained by setting the derivatives of
# E_K[log f_Beta(X)] with respect to alpha and beta to zero.
#
# The resulting equations are:
# psi(alpha) - psi(alpha + beta) = E_K[log(X)]
# psi(beta)  - psi(alpha + beta) = E_K[log(1 - X)].
kumar_to_beta_score <- function(log_ab, omega, dp, rel.tol = 1e-10) {
  
  # Optimize alpha and beta on the log scale to ensure positivity.
  alpha <- exp(log_ab[1])
  beta  <- exp(log_ab[2])
  
  # Expected log terms under Kumaraswamy truth.
  elog_x <- kumar_expected_log_x(omega = omega, dp = dp)
  elog_1mx <- kumar_expected_log1m_x(omega = omega, dp = dp, rel.tol = rel.tol)
  
  # Return the two score equations.
    c(digamma(alpha) - digamma(alpha + beta) - elog_x,
      digamma(beta) - digamma(alpha + beta) - elog_1mx)
}


# ============================================================
# Jacobian of the score equations in log(alpha), log(beta)
# ============================================================

# The Jacobian is computed with respect to log(alpha) and
# log(beta). Since d/dlog(alpha) = alpha * d/dalpha and
# d/dlog(beta) = beta * d/dbeta, the derivatives include
# the corresponding alpha and beta factors.
#
# The derivative of the digamma function is the trigamma
# function, which appears in the Jacobian.
kumar_to_beta_score_jacobian <- function(log_ab) {
  
  alpha <- exp(log_ab[1])
  beta  <- exp(log_ab[2])
  
  tri_sum <- trigamma(alpha + beta)
  
  j11 <- (trigamma(alpha) - tri_sum ) * alpha
  j12 <- -tri_sum * beta
  j21 <- -tri_sum * alpha
  
  j22 <- (trigamma(beta) - tri_sum) * beta
  
  matrix(c(j11, j21,
           j12, j22),
         nrow = 2, ncol = 2)
}


# ============================================================
# Newton solver: Kumaraswamy -> Beta projection
# ============================================================

# Solve the Beta score equations using Newton's method on
# log(alpha) and log(beta).
#
# The Newton update is:
#
# theta_new = theta - J(theta)^(-1) S(theta),
#
# where S is the score vector and J is its Jacobian.
solve_kumar_to_beta <- function(omega, dp, tol = 1e-12, maxit = 100, rel.tol = 1e-10) {
  
  # Obtain starting values from the Kumaraswamy moments.
  log_ab <- kumar_to_beta_start(omega = omega, dp = dp)
  
  converged <- FALSE
  iter <- 0L
  
  for (iter in seq_len(maxit)) {
    
    # Evaluate the score equations at the current parameters.
    score <- kumar_to_beta_score(
      log_ab = log_ab,
      omega = omega,
      dp = dp,
      rel.tol = rel.tol
    )
    
    # Use the Euclidean norm of the score to assess convergence.
    score_norm <- sqrt(sum(score^2))
    
    if (is.finite(score_norm) && score_norm < tol) {
      converged <- TRUE
      break
    }
    
    # Compute the Jacobian of the score equations.
    jac <- kumar_to_beta_score_jacobian(log_ab = log_ab)
    
    # Solve for the Newton step.
    step <- tryCatch(
      solve(jac, score),
      error = function(e) NULL
    )
    
    if (is.null(step) || any(!is.finite(step))) {
      break
    }
    
    # Update log(alpha) and log(beta).
    log_ab_new <- log_ab - step
    
    if (any(!is.finite(log_ab_new)))
      break
    
    log_ab <- log_ab_new
  }
  
  # Recover the positive Beta shape parameters.
  alpha <- unname(exp(log_ab[1]))
  beta <- unname(exp(log_ab[2]))
  
  # Re-evaluate the score at the final solution.
  final_score <- kumar_to_beta_score(
    log_ab = log_ab,
    omega = omega,
    dp = dp, 
    rel.tol = rel.tol
  )
  
  final_score_norm <- sqrt(sum(final_score^2))
  
  if (is.finite(final_score_norm) && final_score_norm < tol) {
    converged <- TRUE
  }
  
  # Convert from shape parameters to the mean-precision
  # parameterization.
  beta_pars <- beta_from_ab(alpha = alpha, beta = beta)
  
  list(
    alpha = alpha,
    beta = beta,
    mu = unname(beta_pars["mu"]),
    phi = unname(beta_pars["phi"]),
    score = unname(final_score),
    score_norm = unname(final_score_norm),
    iterations = as.integer(iter),
    convergence = if (converged) 0L else 1L,
    log_ab = unname(log_ab)
  )
}

# ============================================================
# KL projection: Kumaraswamy -> Beta
# ============================================================

# Compute the Beta KL projection under Kumaraswamy truth.
# This function provides the projection interface to the Newton
# solver and returns the Beta parameters (alpha*, beta*) satisfying
# the expected score equations.

project_kumar_to_beta <- function(omega, dp, rel.tol = 1e-10,
                                  tol = 1e-12, maxit = 100) {
  
  solve_kumar_to_beta(
    omega = omega,
    dp = dp,
    tol = tol,
    maxit = maxit,
    rel.tol = rel.tol
  )
}

# ============================================================
# Expected Kumaraswamy log-density under Kumaraswamy truth
# ============================================================
# 
# Compute the expected log-density of the true Kumaraswamy model:
# E_K[log f_K(X)].
# Both logarithmic expectations have closed-form expressions:
#
# E_K[log(X)] = [psi(1) - psi(q + 1)] / p
# E_K[log(1 - X^p)] = -1 / q.
#
# This quantity is evaluated analytically using the closed-form
# expectations of log(X) and log(1 - X^p).
kumar_expected_log_kumar <- function(omega, dp, rel.tol = 1e-10) {
  
  pars <- kumar_to_pq(omega = omega, dp = dp)
  
  p <- pars["p"]
  q <- pars["q"]
  
  elog_x <- (digamma(1) - digamma(q + 1)) / p
  elog_1m_xp <- -1 / q
  
  log(p) + log(q) + (p - 1) * elog_x + (q - 1) * elog_1m_xp
}


# ============================================================
# Expected Beta log-density under Kumaraswamy truth
# ============================================================
# 
# Compute the expected log-density of a Beta(alpha, beta) model
# under Kumaraswamy truth:
#
# E_K[log f_B(X; alpha, beta)].
#
# When evaluated at the projected parameters (alpha*, beta*),
# this gives the expected log-density of the KL-optimal Beta
# approximation to the Kumaraswamy distribution.
kumar_expected_log_beta <- function(alpha, beta, omega, dp, rel.tol = 1e-10) {
  
  elog_x <- kumar_expected_log_x(omega = omega, dp = dp)
  elog_1mx <- kumar_expected_log1m_x(omega = omega, dp = dp, rel.tol = rel.tol)
  
  -lbeta(alpha, beta) + (alpha - 1) * elog_x + (beta - 1) * elog_1mx
}


# ============================================================
# Full KL projection: Kumaraswamy -> Beta
# ============================================================
# 
# Compute the full KL projection from Kumaraswamy to Beta.
#
# First obtain the KL-optimal Beta parameters (alpha*, beta*)
# from the expected score equations. The intrinsic divergence is
# then evaluated as
#
# D*_(K -> B) = E_K[log f_K(X)] - E_K[log f_B(X; alpha*, beta*)].
#
# The returned D_star therefore measures the intrinsic separation
# from the Kumaraswamy truth to its KL-optimal Beta approximation.
project_kumar_to_beta_kl <- function(omega, dp, rel.tol = 1e-10, 
                                     tol = 1e-12, maxit = 100) {
  
  proj <- project_kumar_to_beta(
    omega = omega,
    dp = dp,
    rel.tol = rel.tol,
    tol = tol,
    maxit = maxit
  )
  
  elog_kumar <- kumar_expected_log_kumar(omega = omega, dp = dp, rel.tol = rel.tol)
  
  elog_beta <- kumar_expected_log_beta(
    alpha = proj$alpha,
    beta = proj$beta,
    omega = omega,
    dp = dp,
    rel.tol = rel.tol
  )
  
  D_star <- unname(elog_kumar - elog_beta)
  
  c(proj,list(D_star = D_star))
}


# ============================================================
# Complete intrinsic separation: Kumaraswamy -> Beta
# ============================================================
#
# Compute the separation between a Kumaraswamy distribution
# and its KL-optimal Beta projection.
#
# The function combines the KL divergence, Hellinger distance,
# overlap coefficient, log-likelihood-ratio variance, and the
# standardized separation delta = D_star / sqrt(v).
# ============================================================
intrinsic_kumar_to_beta <- function(omega, dp, rel.tol = 1e-10, tol = 1e-12,
                                    maxit = 100, zero.tol = 1e-10) {
  
  # Obtain the KL-optimal Beta projection.
  proj <- project_kumar_to_beta_kl(
    omega = omega,
    dp = dp,
    rel.tol = rel.tol,
    tol = tol,
    maxit = maxit
  )
  
  # Define the true Kumaraswamy density and its log-density.
  d_true <- function(x) {
    
    d_kumar(x = x, omega = omega, dp = dp)
  }
  
  log_d_true <- function(x) {
    
    d_kumar(x = x, omega = omega, dp = dp, log = TRUE)
  }
  
  # Define the projected Beta density and its log-density.
  d_comp <- function(x) {
    
    d_beta(x = x, mu = proj$mu, phi = proj$phi)
  }
  
  log_d_comp <- function(x) {
    
    d_beta(x = x, mu = proj$mu, phi = proj$phi, log = TRUE)
  }
  
  # Compute the Hellinger distance between the true and projected densities.
  H_star <- hellinger_distance(log_d_true = log_d_true, log_d_comp = log_d_comp,
                               rel.tol = rel.tol
  )
  
  # Compute the overlap coefficient between the true and projected densities.
  OVL_star <- overlap_coefficient(log_d_true = log_d_true, log_d_comp = log_d_comp,
                                  rel.tol = rel.tol
  )
  
  # KL divergence at the optimal Beta projection.
  D_star <- unname(proj$D_star)
  
  # Variance of the log-likelihood ratio under the Kumaraswamy truth.
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
    truth = "Kumaraswamy",
    omega = unname(omega),
    dp = unname(dp),
    alpha_star = unname(proj$alpha),
    beta_star = unname(proj$beta),
    mu_star = unname(proj$mu),
    phi_star = unname(proj$phi),
    D_star = unname(D_star),
    H_star = unname(H_star),
    OVL_star = unname(OVL_star),
    v = unname(v),
    delta = unname(delta),
    score_norm = unname(proj$score_norm),
    iterations = unname(proj$iterations),
    convergence = unname(proj$convergence)
  )
}
