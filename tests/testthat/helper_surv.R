set.seed(1772)

#' Generate semi-competing risks data and parameter values for local likelihood.
#'
#' @param family Copula family: integer 3-5.
#' @return List with elements `udata`, `epar`, `wgt`, `x`, `x0`, `eta`.
survdata_sim <- function(family) {
  # generate semi-competing risks data
  nobs <- sample(50:100, 1) # number of observations
  x <- sort(runif(nobs)) # values in  (0,1] interval
  # etafun <- function(t) 2*cos(12*pi*t)  # oscillating calibration function
  etafun <- function(t) 2 * t # linear calibration function
  tpar <- LocalCop::BiCopEta2Par(family = family, etafun(x))$par # true copula parameter
  U <- VineCopula::BiCopSim(N = nobs, family = family, par = tpar, par2 = 0)
  mpar1 <- c(1, 1)
  mpar2 <- c(0.5, 0.8)
  mpar3 <- 10
  Y <- -log(U[, 1]) / (mpar1[1] * exp(mpar1[2] * x))
  C <- -log(U[, 2]) / (mpar2[1] * exp(mpar2[2] * x))
  D <- runif(nobs, 0, mpar3)
  L.data <- data.frame(Y, C, D) # true event times
  Y1 <- Y2 <- status1 <- status2 <- status3 <- NULL
  for (i in 1:nobs) {
    Y1[i] <- min(L.data[i, ])
    Y2[i] <- min(L.data[i, -1])
    status1[i] <- ifelse(Y1[i] == Y[i], 1, 0)
    status2[i] <- ifelse(C[i] < D[i], 1, 0)
    status3[i] <- min(1, (status1[i] + status2[i]))
  }
  # get KM estimates of the survival functions (for simplicity in checks)
  u1 <- sapply(Y1, function(y) {
    SurvivalCop::KM(t = y, Y = Y1, status = status1)
  })
  u2 <- sapply(Y2, function(y) {
    SurvivalCop::KM(t = y, Y = Y2, status = status2)
  })
  udata <- data.frame(u1, u2, status1, status2)
  # local likelihood
  x0 <- runif(1, min(x), max(x)) # evaluation point
  # weight specification
  kern <- sample(c(KernEpa, KernGaus, KernBeta, KernBiQuad, KernTriAng), 1)[[1]]
  band <- runif(1, .025, .5)
  wgt <- KernWeight(x = x, x0 = x0, band = band, kernel = kern)
  # local likelihood calculation
  for (ii in 1:100) {
    # generate valid eta/epar pair:
    eta <- rnorm(2) / 2 # evaluation parameter
    epar <- LocalCop::BiCopEta2Par(
      family = family,
      eta = eta[1] + eta[2] * (x - x0)
    )$par
    if (
      (!family %in% c("3", "4")) ||
        ((family == "3") && (max(epar) < 27.9)) ||
        ((family == "4") && (max(epar) < 16.9))
    ) {
      break
    }
  }
  list(udata = udata, epar = epar, wgt = wgt, x = x, x0 = x0, eta = eta)
}


#' Survival copula density
#'
#' @template param-u1
#' @template param-u2
#' @param status1 Vector of censoring indicators for the first variable.
#' @param status2 Vector of censoring indicators for the second variable.
#' @param family Copula family: integer 3-5.
#' @param par Copula parameter value.
#' @return Vector of log-densities.
SurvCopDens <- function(u1, u2, status1, status2, family, par) {
  get.model <- SurvCopModel(family)
  model <- get.model(u1, u2, par)
  n <- length(u1)

  delta1 <- (1 - status1) * (1 - status2)
  delta2 <- status1 * (1 - status2)
  delta3 <- (1 - status1) * status2
  delta4 <- status1 * status2

  res <- sapply(1:n, function(i) {
    term <- ifelse(
      delta1[i] == 1,
      model$cdf[i],
      ifelse(
        delta2[i] == 1,
        model$partial1[i],
        ifelse(delta3[i] == 1, model$partial2[i], model$pdf[i])
      )
    )
    if (is.finite(term)) {
      return(log(term))
    } else {
      return(NA)
    }
  })

  return(res)
}


#' Survival model specification for Archimedean copulas
#'
#' @param family Copula family: integer 3-5.
#' @return List of copula density, h functions and cdf for likelihood calculations.
SurvCopModel <- function(family) {
  if (family == 3) {
    copula_model <- function(u1, u2, par) {
      cdf <- ((u1^(-par)) + (u2^(-par)) - 1)^(-1 / par)
      partial1 <- u1^(-(par + 1)) * (u1^(-par) + u2^(-par) - 1)^(-(1 + 1 / par))
      partial2 <- u2^(-(par + 1)) * (u1^(-par) + u2^(-par) - 1)^(-(1 + 1 / par))
      pdf <- (par + 1) *
        u1^(-(par + 1)) *
        u2^(-(par + 1)) *
        (u1^(-par) + u2^(-par) - 1)^(-(2 + 1 / par))

      model <- list(
        cdf = cdf,
        partial1 = partial1,
        partial2 = partial2,
        pdf = pdf
      )
      return(model)
    }
  }

  if (family == 4) {
    copula_model <- function(u1, u2, par) {
      cdf <- exp(-((-log(u1))^(par) + (-log(u2))^(par))^(1 / (par)))
      partial1 <- (1 / u1) *
        (-log(u1))^((par) - 1) *
        ((-log(u1))^(par) + (-log(u2))^(par))^(1 / (par) - 1) *
        exp(-((-log(u1))^(par) + (-log(u2))^(par))^(1 / (par)))
      partial2 <- (1 / u2) *
        (-log(u2))^((par) - 1) *
        ((-log(u1))^(par) + (-log(u2))^(par))^(1 / (par) - 1) *
        exp(-((-log(u1))^(par) + (-log(u2))^(par))^(1 / (par)))
      pdf <- (1 / u1) *
        (1 / u2) *
        exp(-((-log(u1))^(par) + (-log(u2))^(par))^(1 / (par))) *
        (log(u1) * log(u2))^((par) - 1) *
        ((-log(u1))^(par) + (-log(u2))^(par))^(2 * (1 / (par) - 1)) *
        (1 + ((par) - 1) * ((-log(u1))^(par) + (-log(u2))^(par))^(-1 / (par)))

      model <- list(
        cdf = cdf,
        partial1 = partial1,
        partial2 = partial2,
        pdf = pdf
      )
      return(model)
    }
  }

  if (family == 5) {
    copula_model <- function(u1, u2, par) {
      cdf <- (-1 / par) *
        log(
          1 + (((exp(-par * u1) - 1) * (exp(-par * u2) - 1)) / (exp(-par) - 1))
        )
      partial1 <- exp(-(par) * u1) *
        (exp(-(par) * u2) - 1) /
        ((exp(-(par)) - 1) + (exp(-(par) * u1) - 1) * (exp(-(par) * u2) - 1))
      partial2 <- exp(-(par) * u2) *
        (exp(-(par) * u1) - 1) /
        ((exp(-(par)) - 1) + (exp(-(par) * u1) - 1) * (exp(-(par) * u2) - 1))
      pdf <- ((par) * (1 - exp(-(par))) * exp(-(par) * (u1 + u2))) /
        ((1 - exp(-(par))) - (1 - exp(-(par) * u1)) * (1 - exp(-(par) * u2)))^2

      model <- list(
        cdf = cdf,
        partial1 = partial1,
        partial2 = partial2,
        pdf = pdf
      )
      return(model)
    }
  }
  return(copula_model)
}


#' Survival conditional copula log-likelihood function
#'
#' @template param-u1
#' @template param-u2
#' @param status1 Vector of censoring indicators for the first variable.
#' @param status2 Vector of censoring indicators for the second variable.
#' @param family Copula family: integer 3-5.
#' @param par Copula parameter value.
#' @return Vector of log-densities.
SurvCondiCopLik <- function(
  u1,
  u2,
  status1,
  status2,
  family,
  X,
  eta,
  degree = 1,
  x = NULL,
  band = NULL,
  kernel = KernEpa
) {
  if (!degree %in% 0:1) {
    stop("degree must be 0 or 1.")
  }
  if (!is.null(x) & is.null(band) | is.null(x) & !is.null(band)) {
    stop("Local likelihood calculation requires both x and bandwidth.")
  }

  # define eta
  if (degree == 0) {
    eta_val <- eta[1]
  }
  if (degree == 1) {
    if (!is.null(x)) {
      eta_val <- eta[1] + eta[2] * (X - x)
    } else {
      eta_val <- eta[1] + eta[2] * X
    }
  }
  par_val <- LocalCop::BiCopEta2Par(family = family, eta = eta_val)$par
  if (!is.null(x) & !is.null(band)) {
    wgt <- (1 / band) * kernel((X - x) / band)
    wpos <- wgt > 0
    wgt <- wgt[wpos]
    X <- X[wpos]
    u1 <- u1[wpos]
    u2 <- u2[wpos]
    status1 <- status1[wpos]
    status2 <- status2[wpos]
    par_val <- par_val[wpos]
  }
  logdens <- SurvCopDens(
    u1 = u1,
    u2 = u2,
    status1 = status1,
    status2 = status2,
    family = family,
    par = par_val
  )

  if (!is.null(x) & !is.null(band)) {
    logdens <- logdens * wgt
  }
  loglik <- sum(logdens[is.finite(logdens)])
  return(loglik)
}
