# the following example shows the maximum likelihood estimation
# for constant dependence parameter in bivariate survival data

# simulate survival data
set.seed(51)
n <- 1000 # sample size
family <- 3 # Clayton copula
x <- runif(n)
truetau <- runif(1,0,1)
plot(x, rep(truetau, n), ylim = c(0, 1))
par <- VineCopula::BiCopTau2Par(family, truetau)
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


# (pseudo) maximum likelihood estimation via TMB

system.time({
  eta_hat <- SurvCopFit(
    u1 = u1,
    u2 = u2,
    status1 = status1,
    status2 = status2,
    family = family
  )
})





# optimizing manual likelihood function
fun <- function(v) {
  -SurvCondiCopLik(
    u1 = u1,
    u2 = u2,
    status1 = status1,
    status2 = status2,
    family = family,
    X = 0,
    eta = v, 
    degree = 0)}


system.time({
  eta_hat2 <- nlminb(start=BiCopTau2Eta(family, cor(u1,u2)), objective=fun)$par
})


# compare results
plot(
  x,
  rep(truetau, n),
  type = "l",
  xlab = expression(x),
  ylab = expression(tau(x)),
  ylim = c(0, 1)
)
lines(x, rep(BiCopEta2Tau(family, eta = eta_hat), n), col = "red")
lines(x, rep(BiCopEta2Tau(family, eta = eta_hat2), n), col = "blue")
legend(
  "bottomleft",
  fill = c("black", "red", "blue"),
  legend = c("True", "TMB", "R")
)
