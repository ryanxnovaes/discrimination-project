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