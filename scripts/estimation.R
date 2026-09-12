source("scripts/distributions.R")

# ============================================================
# Beta maximum likelihood estimation
# ============================================================

loglik_beta <- function(mu, phi, x) {
  
  if (any(!is.finite(x)) || any(x <= 0) || any(x >= 1))
    stop("All observations must be finite and lie strictly inside (0, 1)")
  
  sum(d_beta(x = x, mu = mu, phi = phi, log = TRUE))
}

fit_beta <- function(x, start_mu = NULL, start_phi = NULL,
                     method = "L-BFGS-B", maxit = 1000) {
  
  # ----------------------------------------------------------
  # Starting values
  # ----------------------------------------------------------
  
  if (is.null(start_mu))
    start_mu <- mean(x)
  
  if (!is.finite(start_mu) || start_mu <= 0 || start_mu >= 1)
    stop("Invalid starting value for mu")
  
  if (is.null(start_phi)) {
    
    # alpha = mu * phi and beta = (1 - mu) * phi
    #
    # Var(X) = alpha * beta /
    #          ((alpha + beta)^2 * (alpha + beta + 1))
    #
    #        = mu * (1 - mu) / (phi + 1)
    #
    # Therefore:
    #
    # phi = mu * (1 - mu) / Var(X) - 1
    
    phi_mm_value <- start_mu * (1 - start_mu) / stats::var(x) - 1
    
    if (is.finite(phi_mm_value) && phi_mm_value > 0) {
      
      start_phi <- phi_mm_value
      
    } else {
      
      start_phi <- 1
    }
  }
  
  
  # ----------------------------------------------------------
  # Negative log-likelihood
  # ----------------------------------------------------------
  
  objective <- function(par) {
    
    mu <- par["mu"]
    phi <- par["phi"]
    
    -loglik_beta(mu = mu, phi = phi, x = x)
  }
  
  
  # ----------------------------------------------------------
  # Maximum likelihood estimation
  # ----------------------------------------------------------
  
  opt <- optim(
    par = c(mu = start_mu, phi = start_phi),
    fn = objective,
    method = method,
    lower = c(mu = 1e-8, phi = 1e-8),
    upper = c(mu = 1 - 1e-8, phi = Inf),
    control = list(maxit = maxit)
  )
  
  
  # ----------------------------------------------------------
  # Estimates
  # ----------------------------------------------------------
  
  ab_hat <- beta_to_ab(mu = opt$par["mu"], phi = opt$par["phi"])
  
  list(
    mu = opt$par["mu"], 
    phi = opt$par["phi"],
    alpha = ab_hat["alpha"], 
    beta = ab_hat["beta"],
    logLik = -opt$value,
    convergence = opt$convergence,
    counts = opt$counts,
    optim = opt
  )
}


# ============================================================
# Kumaraswamy maximum likelihood estimation
# ============================================================

loglik_kumar <- function(p, q, x) {
  
  if (any(!is.finite(x)) || any(x <= 0) || any(x >= 1))
    stop("All observations must be finite and lie strictly inside (0, 1)")
  
  if (!is.finite(p) || !is.finite(q) || p <= 0 || q <= 0)
    return(-Inf)
  
  log_x <- log(x)
  
  a <- p * log_x
  
  log_one_minus_xp <- log1mexp(a)
  
  ll <- length(x) * (log(p) + log(q)) +
    (p - 1) * sum(log_x) +
    (q - 1) * sum(log_one_minus_xp)
  
  if (!is.finite(ll))
    return(-Inf)
  
  ll
}

kumar_start_pq <- function(x) {
  
  m <- mean(x)
  v <- var(x)
  
  f <- function(par) {
    p <- exp(par[1])
    q <- exp(par[2])
    
    mu <- exp(log(q) + lbeta(1 + 1 / p, q))
    mu2 <- exp(log(q) + lbeta(1 + 2 / p, q))
    
    (mu - m)^2 + (mu2 - mu^2 - v)^2
  }
  
  o <- optim(log(c(1, 1)), 
             f,
             method = "L-BFGS-B",
             lower = log(c(0.01, 0.01)),
             upper = log(c(100, 100)))
  
  c(p = exp(o$par[1]), q = exp(o$par[2]))
}


fit_kumar <- function(x, start_p = NULL, start_q = NULL,
                      method = "BFGS", maxit = 1000) {
  
  if (is.null(start_p) || is.null(start_q)) {
    start <- kumar_start_pq(x)
    if (is.null(start_p)) start_p <- unname(start["p"])
    if (is.null(start_q)) start_q <- unname(start["q"])
  }
  
  
  # ----------------------------------------------------------
  # Negative log-likelihood in unconstrained parameters
  # ----------------------------------------------------------
  
  objective <- function(par) {
    
    p <- exp(par["log_p"])
    q <- exp(par["log_q"])
    
    ll <- loglik_kumar(p = p, q = q, x = x)
    
    -ll
  }
  
  
  # ----------------------------------------------------------
  # Maximum likelihood estimation
  # ----------------------------------------------------------
  
  opt <- optim(
    par = c(log_p = log(start_p), log_q = log(start_q)),
    fn = objective,
    method = method,
    control = list(maxit = maxit))
  
  
  # ----------------------------------------------------------
  # Estimates
  # ----------------------------------------------------------
  
  omega_dp_hat <- kumar_from_pq(p = exp(opt$par["log_p"]), q = exp(opt$par["log_q"]))
  
  list(
    omega = omega_dp_hat["omega"],
    dp = omega_dp_hat["dp"],
    p = unname(exp(opt$par["log_p"])),
    q = unname(exp(opt$par["log_q"])),
    start_p = start_p,
    start_q = start_q,
    logLik = -opt$value,
    convergence = opt$convergence,
    counts = opt$counts,
    optim = opt
  )
}