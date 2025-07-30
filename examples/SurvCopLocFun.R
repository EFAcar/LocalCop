# the following example shows the calculation of an unconditional copula 
# likelihood function for bivariate survival data

# simulate survival data
set.seed(31)
n <- 1000 # sample size
family <- 3 # Clayton copula
rho <- runif(1, 0, 1) # unconditional dependence parameter
par <- VineCopula::BiCopTau2Par(family, rho)
marpar1 <- c(1, 1)
marpar2 <- c(2, 2)
marparc <- c(3, 3)
U <- VineCopula::BiCopSim(N = n, family = family, par = par, par2 = 0)
D <- rweibull(n, shape = marparc[2], scale =  marparc[1]^(-1/marparc[2]))
Y <- qweibull(U[,1], shape = marpar1[2], 
              scale = marpar1[1]^(-1/marpar1[2]), lower.tail= F)
C <- qweibull(U[,2], shape = marpar2[2], 
              scale = marpar2[1]^(-1/marpar2[2]), lower.tail= F)
Y1 <- pmin(Y, D) 
Y2 <- pmin(C, D) 
status1 <- 1*(Y <= D)
status2 <- 1*(C <= D)

# get KM estimates of the survival functions 
km1 <- survival::survfit(survival::Surv(Y1, status1) ~ 1, data=data.frame(Y1, status1))
u1 <- km1$surv[match(Y1, km1$time)]
km2 <- survival::survfit(survival::Surv(Y2, status2) ~ 1, data=data.frame(Y2, status2))
u2 <- km2$surv[match(Y2, km2$time)]

# parameter conversion: equivalent to BiCopPar2Eta(family = 1, ...)
rho2eta <- function(rho) .5 * log((1+rho)/(1-rho))
rhovec <- runif(50, 0, 1)

# create likelihood function for observed data (excluding zero's)
nll_obj <- SurvCopLocFun(u1 = u1, u2 = u2, 
                         status1 = status1, 
                         status2 = status2,
                         family = family, 
                         x = rep(0,n), x0 = 0, # centered covariate x - x0 == 0
                         wgt = rep(1, n),  # unweighted
                         degree = 0, # zero-order fit
                         eta = c(rho2eta(rho), 0), 
                         rm.zero = TRUE)
stucop_lik <- function(rho) {
  -nll_obj$fn(c(rho2eta(rho), 0))
}

system.time({
  ll1 <- sapply(rhovec, stucop_lik) # LocalCop
})



# likelihood function based on SurvivalCop
stucop_lik2 <- function(rho) {
  SurvivalCop::SurvCondiCopLik(u1 = u1, u2 = u2, 
                               status1 = status1, status2 = status2,
                               family = family, 
                               X = 0, eta = c(rho2eta(rho), 0))
}


system.time({
  ll2 <- sapply(rhovec, stucop_lik2) # LocalCop
})

# compare the likelihood functions
plot(ll1, ll2)


