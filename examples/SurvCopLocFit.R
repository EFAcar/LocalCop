# simulate data
set.seed(11)
n <- 1000 # sample size
family <- 3 # Clayton copula
x <- runif(n)
rho <- runif(1, 0, 1) # unconditional dependence parameter
par <- VineCopula::BiCopTau2Par(family, rho)
sdata <- SurvSim(n, family = family, par = par, 
                 marpar1 = c(1, 1), marpar2 = c(2, 2), marparc = c(3, 3),  
                 mardist = "weibull", type = "rcen")

# get observed data
data <- sdata$data
Y1 <- data[,1]
Y2 <- data[,2]
status1 <- data[,3]
status2 <- data[,4]

# get KM estimates of the survival functions 
u1 <- sapply(Y1, function(y){SurvivalCop::KM(t=y, Y = Y1, status = status1)})
u2 <- sapply(Y2, function(y){SurvivalCop::KM(t=y, Y = Y2, status = status2)})


# local likelihood estimation
x0 <- seq(min(x), max(x), len = 100)
band <- .5
system.time({
  eta_hat <- SurvCopLocFit(u1 = u1, u2 = u2,
                           status1 = status1, status2 = status2,
                           family = family, 
                           x = x, x0 = x0, band = band)
})

  
# custom optimization routine using stats::optim (gradient-free)
my_optim <- function(obj) {
  opt <- stats::optim(par = obj$par, fn = obj$fn, method = "Nelder-Mead")
  return(opt$par[1]) # always return constant term, even if degree > 0
}
system.time({
  eta_hat2 <- SurvCopLocFit(u1 = u1, u2 = u2,
                             status1 = status1, status2 = status2,
                             family = family, 
                             x = x, x0 = x0, band = band,
                             optim_fun = my_optim)
})

plot(x0, BiCopEta2Tau(family, eta = eta_fun(x0)), type = "l",
     xlab = expression(x), ylab = expression(tau(x)))
lines(x0, BiCopEta2Tau(family, eta = eta_hat$eta), col = "red")
lines(x0, BiCopEta2Tau(family, eta = eta_hat2$eta), col = "blue")
legend("bottomright", fill = c("black", "red", "blue"),
       legend = c("True", "optim_default", "Nelder-Mead"))
