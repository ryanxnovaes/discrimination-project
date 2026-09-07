# ============================================================
# Beta distribution
# ============================================================

# Convert beta mean-precision parameters to shape parameters
beta_to_ab <- function(mu, phi) {
  
  if (!is.finite(mu) || !is.finite(phi) || mu <= 0 || mu >= 1 || phi <= 0)
    stop("Invalid beta parameters")
  
  mu <- unname(mu)
  phi <- unname(phi)
  
  alpha <- mu * phi
  beta <- (1 - mu) * phi
  
  c(alpha = alpha, beta = beta)
}


# Convert beta shape parameters to mean-precision parameters
beta_from_ab <- function(alpha, beta) {
  
  if (!is.finite(alpha) || !is.finite(beta) || alpha <= 0 || beta <= 0)
    stop("Invalid beta shape parameters")
  
  alpha <- unname(alpha)
  beta <- unname(beta)
  
  phi <- alpha + beta
  mu <- alpha / phi
  
  c(mu = mu, phi = phi)
}


# Beta density reparameterized by mean and precision
d_beta <- function(x, mu, phi, log = FALSE) {
  
  pars <- beta_to_ab(mu, phi)
  
  stats::dbeta(x, shape1 = pars["alpha"], shape2 = pars["beta"], log = log)
}


# Beta CDF reparameterized by mean and precision
p_beta <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  
  pars <- beta_to_ab(mu, phi)
  
  stats::pbeta(q, shape1 = pars["alpha"], shape2 = pars["beta"],
               lower.tail = lower.tail, log.p = log.p)
}


# Beta quantile reparameterized by mean and precision
q_beta <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  
  pars <- beta_to_ab(mu, phi)
  
  stats::qbeta(p, shape1 = pars["alpha"], shape2 = pars["beta"],
               lower.tail = lower.tail, log.p = log.p)
}


# Beta random generator reparameterized by mean and precision
r_beta <- function(n, mu, phi) {
  
  
  if (length(n) != 1L || !is.finite(n) || n <= 0 || n != as.integer(n))
    stop("n must be a positive integer")
  
  pars <- beta_to_ab(mu, phi)
  
  stats::rbeta(n, shape1 = pars["alpha"], shape2 = pars["beta"])
}


# ============================================================
# Kumaraswamy distribution
# ============================================================

if (!requireNamespace("extraDistr", quietly = TRUE)) {
  stop("Package 'extraDistr' must be installed.")
}

# Convert Kumaraswamy omega-dispersion parameters to shape parameters
kumar_to_pq <- function(omega, dp) {
  
  if (!is.finite(omega) || !is.finite(dp) || omega <= 0 || omega >= 1 || dp <= 0)
    stop("Invalid Kumaraswamy parameters")
  
  omega <- unname(omega)
  dp <- unname(dp)
  
  p <- 1 / dp
  q <- log(0.5) / log1p(-omega^(1 / dp))
  
  c(p = p, q = q)
}


# Convert Kumaraswamy shape parameters to omega-dispersion parameters
kumar_from_pq <- function(p, q) {
  
  if (!is.finite(p) || !is.finite(q) || p <= 0 || q <= 0)
    stop("Invalid Kumaraswamy shape parameters")
  
  p <- unname(p)
  q <- unname(q)
  
  dp <- 1 / p
  omega <- (-expm1(log(0.5) / q))^(1 / p)
  
  c(omega = omega, dp = dp)
}

# Kumaraswamy density reparameterized by omega and dispersion
# d_kumar <- function(x, omega, dp, log = FALSE) {
#   
#   pars <- kumar_to_pq(omega = omega, dp = dp)
#   
#   p <- pars["p"]
#   q <- pars["q"]
#   
#   z <- p * log(x)
#   
#   log1mxp <- numeric(length(z))
#   idx <- z < log(0.5)
#   
#   log1mxp[idx] <- log1p(-exp(z[idx]))
#   log1mxp[!idx] <- log(-expm1(z[!idx]))
#   
#   log_density <- log(p) + log(q) +
#     (p - 1) * log(x) +
#     (q - 1) * log1mxp
#   
#   if (log) log_density else exp(log_density)
# }

# Kumaraswamy density reparameterized by omega and dispersion

# ============================================================
# Kumaraswamy distribution
# ============================================================


# Numerically stable evaluation of log(1 - exp(a)), for a <= 0
log1mexp <- function(a) {
  
  if (any(a > 0, na.rm = TRUE))
    stop("All values of 'a' must be less than or equal to zero")
  
  out <- numeric(length(a))
  
  idx <- a < log(0.5)
  
  out[idx] <- log1p(-exp(a[idx]))
  out[!idx] <- log(-expm1(a[!idx]))
  
  out
}


# Convert Kumaraswamy omega-dispersion parameters to shape parameters
kumar_to_pq <- function(omega, dp) {
  
  if (!is.finite(omega) || !is.finite(dp) ||
      omega <= 0 || omega >= 1 || dp <= 0)
    stop("Invalid Kumaraswamy parameters")
  
  omega <- unname(omega)
  dp <- unname(dp)
  
  p <- 1 / dp
  
  q <- log(0.5) /
    log1p(-omega^(1 / dp))
  
  c(p = p, q = q)
}


# Convert Kumaraswamy shape parameters to omega-dispersion parameters
kumar_from_pq <- function(p, q) {
  
  if (!is.finite(p) || !is.finite(q) ||
      p <= 0 || q <= 0)
    stop("Invalid Kumaraswamy shape parameters")
  
  p <- unname(p)
  q <- unname(q)
  
  dp <- 1 / p
  
  omega <- (-expm1(log(0.5) / q))^(1 / p)
  
  c(omega = omega, dp = dp)
}


# ============================================================
# Kumaraswamy density
# ============================================================

d_kumar <- function(x, omega, dp, log = FALSE) {
  
  pars <- kumar_to_pq(
    omega = omega,
    dp = dp
  )
  
  shape_p <- unname(pars["p"])
  shape_q <- unname(pars["q"])
  
  out <- rep(NA_real_, length(x))
  
  
  # ----------------------------------------------------------
  # Outside the support
  # ----------------------------------------------------------
  
  idx_out <- !is.na(x) & (x < 0 | x > 1)
  
  out[idx_out] <- if (log) -Inf else 0
  
  
  # ----------------------------------------------------------
  # Interior of the support
  # ----------------------------------------------------------
  
  idx <- !is.na(x) & x > 0 & x < 1
  
  if (any(idx)) {
    
    log_x <- log(x[idx])
    
    log_one_minus_xp <- log1mexp(
      shape_p * log_x
    )
    
    log_density <-
      log(shape_p) +
      log(shape_q) +
      (shape_p - 1) * log_x +
      (shape_q - 1) * log_one_minus_xp
    
    out[idx] <- if (log)
      log_density
    else
      exp(log_density)
  }
  
  
  # ----------------------------------------------------------
  # Left boundary: x = 0
  # ----------------------------------------------------------
  
  idx_zero <- !is.na(x) & x == 0
  
  if (any(idx_zero)) {
    
    density_zero <- if (shape_p < 1) {
      Inf
    } else if (shape_p == 1) {
      shape_q
    } else {
      0
    }
    
    out[idx_zero] <- if (log)
      log(density_zero)
    else
      density_zero
  }
  
  
  # ----------------------------------------------------------
  # Right boundary: x = 1
  # ----------------------------------------------------------
  
  idx_one <- !is.na(x) & x == 1
  
  if (any(idx_one)) {
    
    density_one <- if (shape_q < 1) {
      Inf
    } else if (shape_q == 1) {
      shape_p
    } else {
      0
    }
    
    out[idx_one] <- if (log)
      log(density_one)
    else
      density_one
  }
  
  out
}


# ============================================================
# Kumaraswamy CDF
# ============================================================

p_kumar <- function(q, omega, dp,
                    lower.tail = TRUE,
                    log.p = FALSE) {
  
  pars <- kumar_to_pq(
    omega = omega,
    dp = dp
  )
  
  shape_p <- unname(pars["p"])
  shape_q <- unname(pars["q"])
  
  out <- rep(NA_real_, length(q))
  
  
  # ----------------------------------------------------------
  # q <= 0
  # ----------------------------------------------------------
  
  idx_lower <- !is.na(q) & q <= 0
  
  if (lower.tail) {
    out[idx_lower] <- if (log.p) -Inf else 0
  } else {
    out[idx_lower] <- if (log.p) 0 else 1
  }
  
  
  # ----------------------------------------------------------
  # q >= 1
  # ----------------------------------------------------------
  
  idx_upper <- !is.na(q) & q >= 1
  
  if (lower.tail) {
    out[idx_upper] <- if (log.p) 0 else 1
  } else {
    out[idx_upper] <- if (log.p) -Inf else 0
  }
  
  
  # ----------------------------------------------------------
  # Interior of the support
  # ----------------------------------------------------------
  
  idx <- !is.na(q) & q > 0 & q < 1
  
  if (any(idx)) {
    
    log_x <- log(q[idx])
    
    # log{1 - x^p}
    log_one_minus_xp <- log1mexp(
      shape_p * log_x
    )
    
    # log survival probability:
    # log{1 - F(x)} = q * log{1 - x^p}
    log_survival <- shape_q * log_one_minus_xp
    
    if (lower.tail) {
      
      if (log.p) {
        
        # log F(x)
        out[idx] <- log1mexp(log_survival)
        
      } else {
        
        # F(x) = 1 - exp(log_survival)
        out[idx] <- -expm1(log_survival)
      }
      
    } else {
      
      if (log.p) {
        
        out[idx] <- log_survival
        
      } else {
        
        out[idx] <- exp(log_survival)
      }
    }
  }
  
  out
}


# ============================================================
# Kumaraswamy quantile function
# ============================================================

q_kumar <- function(p, omega, dp,
                    lower.tail = TRUE,
                    log.p = FALSE) {
  
  pars <- kumar_to_pq(
    omega = omega,
    dp = dp
  )
  
  shape_p <- unname(pars["p"])
  shape_q <- unname(pars["q"])
  
  
  # ----------------------------------------------------------
  # Validate probabilities
  # ----------------------------------------------------------
  
  if (log.p) {
    
    if (any(p > 0, na.rm = TRUE))
      stop("Log-probabilities must be less than or equal to zero")
    
  } else {
    
    if (any((p < 0 | p > 1), na.rm = TRUE))
      stop("Probabilities must lie in [0, 1]")
  }
  
  
  out <- rep(NA_real_, length(p))
  
  idx <- !is.na(p)
  
  if (!any(idx))
    return(out)
  
  
  # ----------------------------------------------------------
  # Compute log(1 - u) without unnecessary cancellation
  #
  # u denotes the lower-tail probability.
  # ----------------------------------------------------------
  
  if (log.p) {
    
    if (lower.tail) {
      
      # p contains log(u)
      log_one_minus_u <- log1mexp(p[idx])
      
    } else {
      
      # p contains log(1 - u)
      log_one_minus_u <- p[idx]
    }
    
  } else {
    
    if (lower.tail) {
      
      # log(1 - u)
      log_one_minus_u <- log1p(-p[idx])
      
    } else {
      
      # p itself is the upper-tail probability = 1 - u
      log_one_minus_u <- log(p[idx])
    }
  }
  
  
  # ----------------------------------------------------------
  # Q(u) = [1 - (1-u)^(1/q)]^(1/p)
  #
  # Stable evaluation:
  #
  # a = log(1-u) / q
  #
  # log{1 - exp(a)} = log1mexp(a)
  # ----------------------------------------------------------
  
  a <- log_one_minus_u / shape_q
  
  log_inner <- log1mexp(a)
  
  out[idx] <- exp(
    log_inner / shape_p
  )
  
  out
}


# ============================================================
# Kumaraswamy random generator
# ============================================================

r_kumar <- function(n, omega, dp) {
  
  if (length(n) != 1L ||
      !is.finite(n) ||
      n <= 0 ||
      n != as.integer(n))
    stop("n must be a positive integer")
  
  u <- stats::runif(n)
  
  q_kumar(
    p = u,
    omega = omega,
    dp = dp
  )
}