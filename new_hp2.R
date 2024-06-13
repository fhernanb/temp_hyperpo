HYPERPO2 <- function (mu.link="log", sigma.link="log") {
  
  mstats <- checklink("mu.link", "HYPERPO2",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "HYPERPO2",
                      substitute(sigma.link), c("log"))
  
  structure(list(family=c("HYPERPO2", "Hyper-Poisson-2"),
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
                 
                 # Primeras derivadas, por ahora son computacionales
                 
                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },
                 
                 # Segundas derivadas, por ahora son computacionales
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },
                 
                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dHYPERPO2(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pHYPERPO2", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),
                 
                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_HYPERPO2(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_HYPERPO2(y)[2], length(y)) ),
                 
                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 
                 y.valid = function(y) all(y >= 0)
                 
  ),
  class=c("gamlss.family", "family"))
}

# obtaining_lambda1 <- function(media, gamma) {
#   # Begin aux function
#   fun <- function(x) x-(gamma-1)*(1-1/f11_cpp(gamma, x))-media
#   fun <- Vectorize(fun)
#   # End aux function
#   if (gamma == 1)
#     result <- media
#   else {
#     res <- uniroot(f=fun,
#                    lower=min(media, max(media+gamma-1, gamma*media)),
#                    upper=max(media, min(media+gamma-1, gamma*media)))
#     result <- res$root
#   }
#   result
# }
# 
# obtaining_lambda2 <- function(media, gamma) {
#   # Begin aux function
#   fun <- function(x) x-(gamma-1)*((F11(gamma, x)-1)/F11(gamma, x))-media
#   fun <- Vectorize(fun)
#   # End aux function
#   if (gamma == 1)
#     result <- media
#   else {
#     res <- uniroot(f=fun,
#                    lower=min(media, max(media+gamma-1, gamma*media)),
#                    upper=max(media, min(media+gamma-1, gamma*media)))
#     result <- res$root
#   }
#   result
# }

obtaining_lambda_cpp <- Vectorize(obtaining_lambda_cpp)

# Chequeo vectorizacion
obtaining_lambda_cpp(media=2, gamma=3)
obtaining_lambda_cpp(media=3, gamma=4)
obtaining_lambda_cpp(media=2:3, gamma=3:4)

dHYPERPO2 <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  
  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda_cpp(media=mu, gamma=sigma)
  
  dHYPERPO(x=x, mu=mu, sigma=sigma, log=log)
}

# Chequeo vectorizacion
dHYPERPO2(x=1, mu=2, sigma=3)
dHYPERPO2(x=2, mu=3, sigma=4)
dHYPERPO2(x=1:2, mu=2:3, sigma=3:4)

pHYPERPO2 <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  
  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda_cpp(media=mu, gamma=sigma)
  
  ly <- max(length(q), length(mu), length(sigma))
  q <- rep(q, length = ly)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)
  
  pHYPERPO(q=q, mu=mu, sigma=sigma, lower.tail=lower.tail, log.p=log.p)
}

# Chequeo vectorizacion
pHYPERPO2(q=0.1, mu=2, sigma=3)
pHYPERPO2(q=0.2, mu=3, sigma=4)
pHYPERPO2(q=c(0.1, 0.2), mu=2:3, sigma=3:4)


rHYPERPO2 <- function(n, mu=1, sigma=1) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))
  
  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda_cpp(media=mu, gamma=sigma)
  
  if (!is.numeric(n) || length(n) != 1 || n < 0)
    stop("invalid arguments")
  if (!(is.double(sigma) || is.integer(sigma)) || !(is.double(mu) ||
                                                    is.integer(mu)))
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
rHYPERPO2(n=1, mu=5, sigma=3)
rHYPERPO2(n=1, mu=50, sigma=4)
rHYPERPO2(n=2, mu=c(2, 30), sigma=3:4)


qHYPERPO2 <- function(p, mu = 1, sigma = 1, lower.tail = TRUE,
                      log.p = FALSE) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda_cpp(media=mu, gamma=sigma)
  
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  qHYPERPO(p=p, mu=mu, sigma=sigma)
}

# Chequeo vectorizacion
qHYPERPO2(p=0.1, mu=5, sigma=3)
qHYPERPO2(p=0.8, mu=6, sigma=4)
qHYPERPO2(p=c(0.1, 0.8), mu=5:6, sigma=3:4)


logLik_HYPERPO2 <- function(logparam=c(0, 0), x){
  return(sum(dHYPERPO2(x     = x,
                       mu    = exp(logparam[1]),
                       sigma = exp(logparam[2]),
                       log=TRUE)))
}

estim_mu_sigma_HYPERPO2 <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_HYPERPO2,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}

