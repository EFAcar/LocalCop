# the following example shows the local likelihood estimation 
# for varying dependence parameter in bivariate survival data

# simulate data
set.seed(31)
n <- 1000 # sample size
family <- 3 # Clayton copula
x <- runif(n)
truetau <- function(x, family){
  eta <- 2*cos(10*x^2-1)-sin(2*x^3)+x
  tau <- BiCopEta2Tau(family = family, eta=eta)
  # tau <- rep(runif(1, 0, 1), length(x)) # constant dependence parameter in Kendall's tau scale
  return(tau)
}
plot(x, truetau(x, family), ylim = c(0,1))
rho <- truetau(x, family)
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

# get KM estimates of the survival functions -- for checks only
u1 <- sapply(Y1, function(y){SurvivalCop::KM(t=y, Y = Y1, status = status1)})
u2 <- sapply(Y2, function(y){SurvivalCop::KM(t=y, Y = Y2, status = status2)})


# local likelihood estimation
x0 <- seq(min(x), max(x), len = 100)
band <- .1
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



# compare results
plot(x0, truetau(x0, family), type = "l", 
     xlab = expression(x), ylab = expression(tau(x)), ylim = c(0,1))
lines(x0, BiCopEta2Tau(family, eta = eta_hat$eta), col = "red")
lines(x0, BiCopEta2Tau(family, eta = eta_hat2$eta), col = "blue")
legend("bottomleft", fill = c("black", "red", "blue"),
       legend = c("True", "optim_default", "Nelder-Mead"))
