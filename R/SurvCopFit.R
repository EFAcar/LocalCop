#' Maximum pseudo-likelihood estimation for survival copulas.
#'
#' Estimate the bivariate copula dependence parameter `eta` for time-to-event outcomes subject to censoring.
#'
#' @template param-u1
#' @template param-u2
#' @param status1 Vector of censoring indicators for the first variable.
#' @param status2 Vector of censoring indicators for the second variable.
#' @template param-family
#' @param eta Optional initial value of the copula dependence parameter (scalar).  If missing will be estimated unconditionally by [VineCopula::BiCopEst()].
#' @param optim_fun Optional specification of likelihood optimization algorithm.  See **Details**.
#' @param rm.zero Logical indicating whether zero values in `u1` and `u2` should be removed when calculating the likelihood.
#' @return
#' \describe{
#'   \item{`eta`}{The estimated dependence parameter under the specified copula family.}
#' }
#' @details By default, optimization is performed with the quasi-Newton algorithm provided by [stats::nlminb()], which uses gradient information provided by automatic differentiation (AD) as implemented by \pkg{TMB}.
#'
#' If the default method is to be overridden, `optim_fun` should be provided as a function taking a single argument corresponding to the output of [SurvCopFun()], and return a scalar value corresponding to the estimate of `eta` at a given covariate value in `x0`.  Note that \pkg{TMB} calculates the *negative* local (log)likelihood, such that the objective function is to be minimized.  See **Examples**.
#' @example examples/SurvCopLocFit.R
#' @export
SurvCopFit <- function(
  u1,
  u2,
  status1,
  status2,
  family,
  eta,
  optim_fun,
  rm.zero = TRUE
) {
  # initialize eta
  etaNu <- .get_etaNu(u1 = u1, u2 = u2, family = family, degree = 0, eta = eta)
  ieta <- etaNu$eta
  # optimization function
  if (missing(optim_fun)) {
    optim_fun <- .optim_default
  }

  wgt <- rep(1, length(u1))
  obj <- SurvCopLocFun(
    u1 = u1,
    u2 = u2,
    status1 = status1,
    status2 = status2,
    family = family,
    x = wgt,
    x0 = 0,
    wgt = wgt,
    degree = 0,
    eta = ieta,
    rm.zero = rm.zero
  )
  eta0 <- optim_fun(obj)
  return(eta = as.numeric(eta0))
}
