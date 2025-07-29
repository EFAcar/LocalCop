#' Local likelihood estimation for survival copulas.
#'
#' Estimate the bivariate copula dependence parameter `eta` at multiple covariate values for time-to-event outcomes subject to censoring. 
#'
#' @template param-u1
#' @template param-u2
#' @param status1 Vector of censoring indicators for the first variable.
#' @param status2 Vector of censoring indicators for the second variable.
#' @template param-family
#' @template param-x
#' @template param-xseq
#' @param nx If `x0` is missing, defaults to `nx` equally spaced values in `range(x)`.
#' @template param-degree
#' @param eta Optional initial value of the copula dependence parameter (scalar).  If missing will be estimated unconditionally by [VineCopula::BiCopEst()].
#' @template param-kernel
#' @template param-band
#' @param optim_fun Optional specification of local likelihood optimization algorithm.  See **Details**.
#' @param rm.zero Logical indicating whether zero values in `u1` and `u2` should be removed when calculating the local likelihood. 
#' @param cl Optional parallel cluster created with [parallel::makeCluster()], in which case optimization for each element of `x0` will be done in parallel on separate cores.  If `cl == NA`, computations are run serially.
#' @return List with the following elements:
#' \describe{
#'   \item{`x`}{The vector of covariate values `x0` at which the local likelihood is fit.}
#'   \item{`eta`}{The vector of estimated dependence parameters of the same length as `x0`.}
#' }
#' @details By default, optimization is performed with the quasi-Newton algorithm provided by [stats::nlminb()], which uses gradient information provided by automatic differentiation (AD) as implemented by \pkg{TMB}.
#'
#' If the default method is to be overridden, `optim_fun` should be provided as a function taking a single argument corresponding to the output of [SurvCopLocFun()], and return a scalar value corresponding to the estimate of `eta` at a given covariate value in `x0`.  Note that \pkg{TMB} calculates the *negative* local (log)likelihood, such that the objective function is to be minimized.  See **Examples**.
#' @example examples/SurvCopLocFit.R
#' @export
SurvCopLocFit <- function(u1, u2, status1, status2, 
                          family, x, x0, nx = 100,
                          degree = 1,
                          eta, kernel = KernEpa, band,
                          optim_fun, rm.zero = TRUE, cl = NA) {

  # default x0
  if(missing(x0)) {
    x0 <- seq(min(x), max(x), len = nx)
  } else {
    x0 <- sort(x0)
  }
  nx <- length(x0)
  # initialize eta and nu
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  ## degree <- match.arg(degree)
  etaNu <- .get_etaNu(u1 = u1, u2 = u2, family = family,
                      degree = degree, eta = eta)
  ieta <- etaNu$eta
  # optimization function
  if(missing(optim_fun)) {
    optim_fun <- .optim_default
  }
  fun <- function(xi) {
    wgt <- KernWeight(x = x, x0 = xi, band = band,
                      kernel = kernel, band_type = "constant")
    obj <- SurvCopLocFun(u1 = u1, u2 = u2, 
                         status1 = status1, status2 = status2,
                         family = family, x = x, x0 = xi,
                         wgt = wgt, degree = degree, eta = ieta, 
                         rm.zero = rm.zero)
    return(optim_fun(obj))
  }
  if(nx == 1) {
    eta0 <- fun(x0)
  } else {
    if(!.check_parallel(cl)) {
      # run serially
      eta0 <- sapply(x0, fun)
    } else {
      # run in parallel
      parallel::clusterExport(cl,
                              varlist = c("fun", "u1", "u2", "status1", "status2",
                                          "family", "x", "band", "kernel", 
                                          "optim_fun","ieta", "rm.zero"),
                              envir = environment())
      eta0 <- parallel::parSapply(cl, X = x0, FUN = fun)
    }
  }
  return(list(x = x0, eta = as.numeric(eta0)))
}

