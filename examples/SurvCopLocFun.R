# the following example shows the calculation of 
# an unconditional copula likelihood function for survival data

# simulate data
set.seed(11)
n <- 1000 # sample size
family <- 3 # Clayton copula
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

# add a small quantity to 'zero' survival estimates
u1_new <- ifelse(u1==0, u1 + runif(1, min = 0, max = 10^-5), u1)
u2_new <- ifelse(u2==0, u2 + runif(1, min = 0, max = 10^-5), u2)

# parameter conversion: equivalent to BiCopPar2Eta(family = 1, ...)
rho2eta <- function(rho) .5 * log((1+rho)/(1-rho))
rhovec <- runif(50, 0, 1)


# create likelihood function for observed data (excluding zero's)
ix <- which(u1==0 | u2==0)
nll_obj <- SurvCopLocFun(u1 = u1, u2 = u2, 
                         status1 = status1, 
                         status2 = status2,
                         family = family, 
                         x = rep(0,n), x0 = 0, # centered covariate x - x0 == 0
                         wgt = rep(1, n), # unweighted
                         degree = 0, # zero-order fit
                         eta = c(rho2eta(rho), 0), 
                         rm.zero = FALSE)
stucop_lik <- function(rho) {
  -nll_obj$fn(c(rho2eta(rho), 0))
}

system.time({
  ll1 <- sapply(rhovec, stucop_lik) # LocalCop
})

# create likelihood function for observed data (with slightly elevated zero's)
nll_obj_new <- SurvCopLocFun(u1 = u1_new, u2 = u2_new, 
                         status1 = status1, status2 = status2,
                         family = family, 
                         x = rep(0,n), x0 = 0, # centered covariate x - x0 == 0
                         wgt = rep(1, n), # unweighted
                         degree = 0, # zero-order fit
                         eta = c(rho2eta(rho), 0), 
                         rm.zero = FALSE)
stucop_lik_new <- function(rho) {
  -nll_obj_new$fn(c(rho2eta(rho), 0))
}

system.time({
  ll1_new <- sapply(rhovec, stucop_lik_new) # LocalCop
})



# likelihood function based on SurvivalCop
stucop_lik2 <- function(rho) {
  SurvivalCop::SurvCondiCopLik(u1 = u1, u2 = u2, 
                               status1 = status1, status2 = status2,
                               family = family, 
                               X = 0, eta = c(rho2eta(rho), 0))
}

stucop_lik2_new <- function(rho) {
  SurvivalCop::SurvCondiCopLik(u1 = u1_new, u2 = u2_new, 
                               status1 = status1, status2 = status2,
                               family = family, 
                               X = 0, eta = c(rho2eta(rho), 0))
} 

system.time({
  ll2 <- sapply(rhovec, stucop_lik2) # LocalCop
})

system.time({
ll2_new <- sapply(rhovec, stucop_lik2_new) # LocalCop
})

# compare the estimates
range(ll1 - ll2)
cbind(ll1,ll2, ll1_new, ll2_new, ll1-ll1_new, ll2-ll2_new)
plot(ll1, ll2) 
points(ll1_new, ll2_new, col=2)

