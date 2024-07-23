#--- test survival likelihood implementation in TMB -------------------------------

## library(LocalCop)
## library(TMB)
## library(testthat)
## library(SurvivalCop)
## source("tests/testthat/helper_surv.R") 

test_that("Survival likelihood is same in manual calculation and TMB", {
  nreps <- 20
  test_descr <- expand.grid(family = 3, # only Archimedean families
                            stringsAsFactors = FALSE)
  n_test <- nrow(test_descr)
  for(ii in 1:n_test) {
    for(jj in 1:nreps) {
      # generate data
      family <- test_descr$family[ii]
      args <- survdata_sim(family = family)
      # args <- data_checks(args)
      # gather data
      u1 <- args$udata[,1]
      u2 <- args$udata[,2]
      status1 <- args$udata[,3]
      status2 <- args$udata[,4]
      epar <- args$epar
      x <- args$x
      wgt <- args$wgt
      wpos <- wgt > 0 # index of positive weights
      # subset data
      wgt <- wgt[wpos]
      x <- x[wpos]
      u1 <- u1[wpos]
      u2 <- u2[wpos]
      status1 <- status1[wpos]
      status2 <- status2[wpos]
      epar <- epar[wpos]
      # censoring groups
      delta1 <- (1-status1)*(1-status2) # both censored
      delta2 <- status1*(1-status2) # first uncensored, second censored
      delta3 <- (1-status1)*status2 # first censored, second uncensored
      delta4 <- status1*status2 # both uncensored
      ix <- c(which(delta4==1), which(delta3==1), which(delta2==1), which(delta1==1))
      # reorder data based on censoring groups
      wgt <- wgt[ix]
      x <- x[ix]
      u1 <- u1[ix]
      u2 <- u2[ix]
      status1 <- status1[ix]
      status2 <- status2[ix]
      epar <- epar[ix]
      # loglik in R
      ll_r <- SurvCopDens1(u1 = u1, 
                           u2 = u2,
                           status1 = status1, 
                           status2 = status2, 
                           family = family,
                           par = epar)
      ll_r <- sum(wgt * ll_r)
      # loglik in TMB
      ll_tmb <- SurvCopLocFun(u1 = u1, 
                              u2 = u2,
                              status1 = status1, 
                              status2 = status2, 
                              family = family, 
                              x = x, 
                              x0 = args$x0,
                              wgt = wgt, 
                              eta = args$eta)
      ll_tmb <- -ll_tmb$fn(args$eta)
      expect_equal(ll_r, ll_tmb)
    }
  }
})
