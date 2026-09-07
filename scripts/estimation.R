source("scripts//distributions.R")

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
    
    start_phi <- start_mu * (1 - start_mu) / stats::var(x) - 1
  }
  
  if (!is.finite(start_phi) || start_phi <= 0)
    stop("Invalid method-of-moments starting value for phi")
  
  
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