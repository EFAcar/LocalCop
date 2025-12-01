# the following example shows the local likelihood estimation
# for varying dependence parameter in bivariate survival data

# simulate survival data
set.seed(51)
n <- 200 # sample size
family <- 3 # Clayton copula
x <- runif(n)
truetau <- function(x, family) {
  eta <- 2 * cos(10 * x^2 - 1) - sin(2 * x^3) + x
  tau <- BiCopEta2Tau(family = family, eta = eta)
  # tau <- rep(runif(1, 0, 1), length(x)) # constant dependence parameter in Kendall's tau scale
  return(tau)
}
plot(x, truetau(x, family), ylim = c(0, 1))
rho <- truetau(x, family)
par <- VineCopula::BiCopTau2Par(family, rho)
marpar1 <- c(1, 1)
marpar2 <- c(2, 2)
marparc <- c(3, 3)
U <- VineCopula::BiCopSim(N = n, family = family, par = par, par2 = 0)
D <- rweibull(n, shape = marparc[2], scale = marparc[1]^(-1 / marparc[2]))
Y <- qweibull(
  U[, 1],
  shape = marpar1[2],
  scale = marpar1[1]^(-1 / marpar1[2]),
  lower.tail = F
)
C <- qweibull(
  U[, 2],
  shape = marpar2[2],
  scale = marpar2[1]^(-1 / marpar2[2]),
  lower.tail = F
)
Y1 <- pmin(Y, D)
Y2 <- pmin(C, D)
status1 <- 1 * (Y <= D)
status2 <- 1 * (C <= D)

# get KM estimates of the survival functions
km1 <- survival::survfit(
  survival::Surv(Y1, status1) ~ 1,
  data = data.frame(Y1, status1)
)
u1 <- km1$surv[match(Y1, km1$time)]
km2 <- survival::survfit(
  survival::Surv(Y2, status2) ~ 1,
  data = data.frame(Y2, status2)
)
u2 <- km2$surv[match(Y2, km2$time)]


# cv likelihood calculation 
band <- seq(0.1, 0.5, by=0.025)
system.time({
  cv_lik <- sapply(band, function(h){
  LocalCop::SurvCopLikCV(u1, u2, status1, status2, family, x,  band = h)})
  })
plot(band, cv_lik)

# local likelihood estimation
x0 <- seq(min(x), max(x), len = 100)
band <- band[which.max(cv_lik)]
system.time({
  eta_hat <- SurvCopLocFit(
    u1 = u1,
    u2 = u2,
    status1 = status1,
    status2 = status2,
    family = family,
    x = x,
    x0 = x0,
    band = band
  )
})


# custom optimization routine using stats::optim (gradient-free)
my_optim <- function(obj) {
  opt <- stats::optim(par = obj$par, fn = obj$fn, method = "Nelder-Mead")
  return(opt$par[1]) # always return constant term, even if degree > 0
}
system.time({
  eta_hat2 <- SurvCopLocFit(
    u1 = u1,
    u2 = u2,
    status1 = status1,
    status2 = status2,
    family = family,
    x = x,
    x0 = x0,
    band = band,
    optim_fun = my_optim
  )
})


# compare results
plot(
  x0,
  truetau(x0, family),
  type = "l",
  xlab = expression(x),
  ylab = expression(tau(x)),
  ylim = c(0, 1)
)
lines(x0, BiCopEta2Tau(family, eta = eta_hat$eta), col = "red")
lines(x0, BiCopEta2Tau(family, eta = eta_hat2$eta), col = "blue")
legend(
  "bottomleft",
  fill = c("black", "red", "blue"),
  legend = c("True", "optim_default", "Nelder-Mead")
)
