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
  
  mu <- alpha / phi
  phi <- alpha + beta
  
  
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
d_kumar <- function(x, omega, dp, log = FALSE) {
  
  pars <- kumar_to_pq(omega = omega, dp = dp)
  
  extraDistr::dkumar(x, a = pars["p"], b = pars["q"], log = log)
}


# Kumaraswamy CDF reparameterized by omega and dispersion
p_kumar <- function(q, omega, dp, lower.tail = TRUE, log.p = FALSE) {
  
  pars <- kumar_to_pq(omega = omega, dp = dp)
  
  extraDistr::pkumar(q, a = pars["p"], b = pars["q"],
                     lower.tail = lower.tail, log.p = log.p)
}


# Kumaraswamy quantile reparameterized by omega and dispersion
q_kumar <- function(p, omega, dp, lower.tail = TRUE, log.p = FALSE) {
  
  pars <- kumar_to_pq(omega = omega, dp = dp)
  
  extraDistr::qkumar(p, a = pars["p"], b = pars["q"],
    lower.tail = lower.tail, log.p = log.p)
}


# Kumaraswamy random generator reparameterized by omega and dispersion
r_kumar <- function(n, omega, dp) {
  
  if (length(n) != 1L || !is.finite(n) || n <= 0 || n != as.integer(n))
    stop("n must be a positive integer")
  
  pars <- kumar_to_pq(omega = omega, dp = dp)
  
  extraDistr::rkumar(n, a = pars["p"], b = pars["q"])
}