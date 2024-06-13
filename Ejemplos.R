# Loading the gamlss and functions
library(gamlss)

source("new_hp.R")
source("new_hp2.R")

# -------------------------------------------------------------------------
# ---------------------------- Part 1 -------------------------------------
# -------------------------------------------------------------------------


# Estimating mu and sigma for HYPERPO -------------------------------------

true_mu <- 1.75
true_sigma <- 2

y <- rHYPERPO(n=5000, mu=true_mu, sigma=true_sigma)

# Using gamlss
mod1 <- NULL
mod1 <- gamlss(y~1, family=HYPERPO,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# The estimated coefficients
c(exp(coef(mod1, what="mu")), exp(coef(mod1, what="sigma")))


# Estimating parameters for a regression model ----------------------------

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(1.21 - 3 * x1) # 0.75 en promedio
  sigma <- exp(1.26 - 2 * x2) # 1.30 en promedio
  y <- rHYPERPO(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2, mu=mu, sigma=sigma)
}

datos <- gendat(n=400)

mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=HYPERPO, data=datos,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)

# -------------------------------------------------------------------------
# ---------------------------- Part 2 -------------------------------------
# -------------------------------------------------------------------------


# Estimating mu and sigma for HYPERPO2 ------------------------------------

true_mu <- 1.75
true_sigma <- 2

y <- rHYPERPO2(n=500, mu=true_mu, sigma=true_sigma)

# Using gamlss
mod3 <- NULL
mod3 <- gamlss(y~1, family=HYPERPO2,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# The estimated coefficients
c(exp(coef(mod3, what="mu")), exp(coef(mod3, what="sigma")))


# Estimating parameters for a regression model ----------------------------

library(gamlss)

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(1.21 - 3 * x1) # 0.75 en promedio
  sigma <- exp(1.26 - 2 * x2) # 1.30 en promedio
  y <- rHYPERPO2(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2, mu=mu, sigma=sigma)
}

datos <- gendat(n=100)

mod4 <- NULL
mod4 <- gamlss(y~x1, sigma.fo=~x2, family=HYPERPO2, data=datos,
               control=gamlss.control(n.cyc=500, trace=TRUE))

library(DGLMExtPois)
fit <- glm.hP(formula.mu = y ~ x1,
              formula.gamma = y ~ x2, 
              data = datos)

# Comparing our proposal versus DGLMEPois
summary(mod4)
summary(fit)

coefs_our <- c(coef(mod4, "mu"), coef(mod4, "sigma"))
coefs_competitor <- unlist(coef(fit))
cbind(true=c(1.21, -3, 1.26, -2), coefs_our, coefs_competitor)


# Comparing processing times between our proposal and DGLMEPois -----------


library(microbenchmark)

microbenchmark(
{
  datos <- gendat(n=80)
  gamlss(y~x1, sigma.fo=~x2, family=HYPERPO2, data=datos,
  control=gamlss.control(n.cyc=500, trace=FALSE))
}, 
{
  datos <- gendat(n=80)
  glm.hP(formula.mu = y ~ x1,
         formula.gamma = y ~ x2, 
         data = datos)
  },
times=10)


