set.seed(1772)

#' Generate semi-competing risks data and parameter values for local likelihood.
#'
#' @param family Copula family: integer 3-5.
#' @return List with elements `udata`, `epar`, `wgt`, `x`, `x0`, `eta`.
survdata_sim <- function(family) {
  # generate semi-competing risks data
  nobs <- sample(50:100, 1)  # number of observations
  x <- sort(runif(nobs))    # values in  (0,1] interval
  # etafun <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
  etafun <- function(t) 2*t   # linear calibration function
  tpar <- LocalCop::BiCopEta2Par(family = family, etafun(x))$par # true copula parameter
  U <- VineCopula::BiCopSim(N = nobs, family = family, par = tpar, par2 = 0)
  mpar1 = c(1, 1) 
  mpar2 = c(0.5, 0.8)
  mpar3 = 10
  Y <- -log(U[, 1])/(mpar1[1]*exp(mpar1[2]*x))
  C <- -log(U[, 2])/(mpar2[1]*exp(mpar2[2]*x))
  D <-  runif(nobs, 0, mpar3)
  L.data <- data.frame(Y, C, D)  # true event times
  Y1 <- Y2 <- status1 <- status2 <- status3 <- NULL
  for(i in 1:nobs){
    Y1[i] <- min(L.data[i, ])
    Y2[i] <- min(L.data[i, -1])
    status1[i] <- ifelse(Y1[i] == Y[i], 1, 0)
    status2[i] <- ifelse(C[i] < D[i], 1, 0)
    status3[i] <- min(1, (status1[i] + status2[i]))
  }
  # get KM estimates of the survival functions (for simplicity in checks)
  u1 <- sapply(Y1, function(y){SurvivalCop::KM(t=y, Y = Y1, status = status1)})
  u2 <- sapply(Y2, function(y){SurvivalCop::KM(t=y, Y = Y2, status = status2)})
  udata <- data.frame(u1, u2, status1, status2)  
  # local likelihood 
  x0 <- runif(1, min(x), max(x)) # evaluation point
  # weight specification
  kern <- sample(c(KernEpa, KernGaus, KernBeta,
                   KernBiQuad, KernTriAng), 1)[[1]]
  band <- runif(1, .025, .5)
  wgt <- KernWeight(x = x, x0 = x0, band = band, kernel = kern)
  # local likelihood calculation
  for(ii in 1:100) {
    # generate valid eta/epar pair:
    eta <- rnorm(2)/2  # evaluation parameter
    epar <- LocalCop::BiCopEta2Par(family = family, eta = eta[1] + eta[2] * (x-x0))$par
    if((!family %in% c("3", "4")) ||
       ((family == "3") && (max(epar) < 27.9)) ||
       ((family == "4") && (max(epar) < 16.9))) {
      break
    }
  }
  list(udata = udata,
       epar = epar, wgt = wgt,
       x = x, x0 = x0, eta = eta)
}


data_checks <- function(args){
  # keep only positive weights
  udata <- args$udata[args$wgt>0,]
  epar <- args$epar[args$wgt>0]
  wgt <- args$wgt[args$wgt>0]
  x <- args$x[args$wgt>0]
  # only keep points with positive survival estimates (exclude zero)
  ix <- udata[,1] > 0 | udata[,2] > 0
  udata <- udata[ix,]
  epar <- epar[ix]
  wgt <- wgt[ix]
  x <- x[ix]
  # return list
  list(udata = udata,
       epar = epar, wgt = wgt,
       x = x, x0 = args$x0, eta = args$eta)
}
