HYPERPO <- function (mu.link="log", sigma.link="log") {
  
  mstats <- checklink("mu.link", "HYPERPO",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "HYPERPO",
                      substitute(sigma.link), c("log"))
  
  structure(list(family=c("HYPERPO", "Hyper-Poisson"),
                 parameters=list(mu=TRUE, sigma=TRUE),
                 nopar=2,
                 type="Discrete",
                 
                 mu.link    = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 
                 mu.linkfun    = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 
                 mu.linkinv    = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr    = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 
                 # First derivates
                 
                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },
                 
                 # Second derivates
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },
                 
                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dHYPERPO(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pHYPERPO", type="Discrete",
                                               ymin=0, y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_HYPERPO(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_HYPERPO(y)[2], length(y)) ),
                 
                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 
                 y.valid = function(y) all(y >= 0)
                 
  ),
  class=c("gamlss.family", "family"))
}

library(Rcpp)
sourceCpp("F11.cpp")

# dHYPERPO <- function(x, mu=1, sigma=1, log=FALSE){
#   if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
#   if (any(mu <= 0))     stop("parameter mu has to be positive!")
#   
#   p1 <- x * log(mu) - lgamma(sigma+x) + lgamma(sigma)
#   f11 <- f11_cpp(gamma=sigma, lambda=mu)
#   p2 <- log(f11)
#   res <- p1 - p2
#   res[x < 0] <- -Inf
#   if(log)
#     return(res)
#   else
#     return(exp(res))
# }

#dHYPERPO <- Vectorize(dHYPERPO)

dHYPERPO <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  
  temp <- cbind(x, mu, sigma, log)
  dHYPERPO_vec(x=temp[, 1], mu=temp[, 2], sigma=temp[, 3], log=temp[,4])
}

# Chequeo vectorizacion
dHYPERPO(x=1, mu=2, sigma=3)
dHYPERPO(x=2, mu=3, sigma=4)
dHYPERPO(x=1:2, mu=2:3, sigma=3:4)

pHYPERPO <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  
  ly <- max(length(q), length(mu), length(sigma))
  q <- rep(q, length = ly)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)
  # Begin auxiliar function
  aux_func <- function(q, mu, sigma) {
    cdf <- numeric(length(q))
    for (i in 1:length(q)) {
      res <- dHYPERPO(x=-1:q[i], mu=mu[i], sigma=sigma[i], log=FALSE)
      cdf[i] <- sum(res)
    }
    cdf
  }
  # End auxiliar function
  cdf <- aux_func(q=q, mu=mu, sigma=sigma)
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

# Chequeo vectorizacion
pHYPERPO(q=0.1, mu=2, sigma=3)
pHYPERPO(q=0.2, mu=3, sigma=4)
pHYPERPO(q=c(0.1, 0.2), mu=2:3, sigma=3:4)


simulate_hp <- function(sigma, mu) {
  pochammer <- function(a, r) if (r == 0) 1 else prod(a:(a + r - 1))
  u <- stats::runif(1)
  y <- 0
  p <- 0
  value <- f11_cpp(gamma=sigma, lambda=mu)
  while (p < u) {
    p <- p + mu ^ y / (value * pochammer(sigma, y))
    y <- y + 1
  }
  y - 1
}

rHYPERPO <- function(n, mu=1, sigma=1) {
  if (!is.numeric(n) || length(n) != 1 || n < 0) 
    stop("invalid arguments")
  if (!(is.double(sigma) || is.integer(sigma)) || !(is.double(mu) || is.integer(mu))) 
    stop("Non-numeric argument to mathematical function")
  sigma <- rep(sigma, length.out = n)
  mu <- rep(mu, length.out = n)
  result <- numeric(length = n)
  warn <- FALSE
  for (ind in seq_len(n)) {
    if (sigma[ind] <= 0 || mu[ind] <= 0) {
      result[ind] <- NaN
      warn <- TRUE
    }
    else {
      result[ind] <- simulate_hp(sigma[ind], mu[ind])
    }
  }
  if (warn) 
    warning("NaN(s) produced: sigma and mu must be strictly positive")
  result
}

# Chequeo vectorizacion
rHYPERPO(n=1, mu=5, sigma=3)
rHYPERPO(n=1, mu=50, sigma=4)
rHYPERPO(n=2, mu=c(2, 30), sigma=3:4)

qHYPERPO <- function(p, mu = 1, sigma = 1, lower.tail = TRUE,
                     log.p = FALSE) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  # Begin auxiliar function
  one_quantile_hyperpo <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      prob <- dHYPERPO(x=0, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dHYPERPO(x=i, mu=mu, sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_hyperpo <- Vectorize(one_quantile_hyperpo)
  # End auxiliar function
  one_quantile_hyperpo(p=p, mu=mu, sigma=sigma)
}

# Chequeo vectorizacion
qHYPERPO(p=0.1, mu=5, sigma=3)
qHYPERPO(p=0.8, mu=6, sigma=4)
qHYPERPO(p=c(0.1, 0.8), mu=5:6, sigma=3:4)


F11 <- function(gamma, lambda, maxiter_series = 10000, tol = 1.0e-10) {
  fac  <- 1
  temp <- 1
  L    <- gamma
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / L
    series <- temp + fac
    if (stopping(series - temp, tol)){
      return(Re(series))
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
  return(Re(series))
}

stopping <- function (x, tol) {
  all(abs(x) <= tol, na.rm = TRUE)
}


logLik_HYPERPO <- function(logparam=c(0, 0), x){
  return(sum(dHYPERPO(x     = x,
                      mu    = exp(logparam[1]),
                      sigma = exp(logparam[2]),
                      log=TRUE)))
}
#' Initial values for hyper Poisson
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @keywords internal
#' @export
#' @importFrom stats optim
estim_mu_sigma_HYPERPO <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_HYPERPO,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}

