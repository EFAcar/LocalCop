# the following example shows the calculation of 
# an unconditional copula likelihood function for survival data

# simulate data
n <- 1000 # sample size
family <- 3 # Clayton copula
rho <- runif(1, -1, 1) # unconditional dependence parameter
par <- VineCopula::BiCopTau2Par(family, rho)
sdata <- SurvSim(n, family = family, par = par, 
                 marpar1 = c(1, 1), marpar2 = c(2, 2), marparc = c(3, 3),  
                 mardist = "weibull", type = "rcen")

# create likelihood function

# parameter conversion: equivalent to BiCopPar2Eta(family = 2, ...)
rho2eta <- function(rho) .5 * log((1+rho)/(1-rho))
nll_obj <- CondiCopLocFun(u1 = udata[,1], u2 = udata[,2], family = family,
                          x = rep(0, n), x0 = 0, # centered covariate x - x0 == 0
                          wgt = rep(1, n), # unweighted
                          degree = 0, # zero-order fit
                          eta = c(rho2eta(rho), 0),
                          nu = nu)

# likelihood function: recall that TMB requires a _negative_ ll
stucop_lik <- function(rho) {
  -nll_obj$fn(c(rho2eta(rho), 0))
}

# compare to VineCopula.
rhovec <- runif(50, -1, 1)
system.time({
  ll1 <- sapply(rhovec, stucop_lik) # LocalCop
})
system.time({
  ll2 <- sapply(rhovec, function(rho) {
    # VineCopula
    sum(log(VineCopula::BiCopPDF(u1 = udata[,1], u2 = udata[,2],
                                 family = family,
                                 par = rho, par2 = nu)))
  })
})

# difference between the two
range(ll1 - ll2)
