#' Cross-validated survival copula likelihood.
#'
#' Leave-one-out local likelihood copula parameter estimates are obtained and used to calculate the conditional copula likelihood function under censoring. 
#'
#' @template param-u1
#' @template param-u2
#' @param status1 Vector of censoring indicators for the first variable.
#' @param status2 Vector of censoring indicators for the second variable.
#' @template param-family
#' @template param-x
#' @template param-degree
#' @param eta,kernel,band,optim_fun,cl See [SurvCopLocFit()].
#' @template param-cv_all
#' @param cveta_out If `TRUE`, return the CV estimate of eta at each point in `x` in addition to the CV log-likelihood.
#' @return If `cveta_out = FALSE`, scalar value of the cross-validated log-likelihood.  Otherwise, a list with elements:
#' \describe{
#'   \item{`x`}{The sorted values of `x`.}
#'   \item{`eta`}{The leave-one-out estimates interpolated from the values in `xind` to all of those in `x`.}
#'   \item{`loglik`}{The cross-validated log-likelihood.}
#' }
#' @export
SurvCopLikCV <- function(u1, u2, status1, status2, family, x, 
                         degree = 1, eta, kernel = KernEpa, band,
                         optim_fun, cveta_out = FALSE,
                         cv_all = FALSE, cl = NA) {

  # sort observations
  ix <- order(x)
  x <- x[ix]
  u1 <- u1[ix]
  u2 <- u2[ix]
  status1 <- status1[ix]
  status2 <- status2[ix]
  
  # initialize eta
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  ## degree <- match.arg(degree)
  ieta <- LocalCop:::.get_etaNu(u1 = u1, u2 = u2, family = family,
                                degree = degree, eta = eta)$eta
  
  # cross validation: estimation step
  # optimization function
  if(missing(optim_fun)) {
    optim_fun <- .optim_default
  }
  fun <- function(ii) {
    wgt <- KernWeight(x = x[-ii], x0 = x[ii], band = band,
                      kernel = kernel, band_type = "constant")
    obj <- SurvCopLocFun(u1 = u1[-ii], u2 = u2[-ii], 
                         status1 = status1[-ii], status2 = status2[-ii], 
                         family = family, x = x[-ii], x0 = x[ii],
                         wgt = wgt, degree = degree, eta = ieta)
    return(optim_fun(obj))
  }
  if(!.check_parallel(cl)) {
    # run serially
    cveta <- sapply(seq_along(x), fun)
  } else {
    # run in parallel
    parallel::clusterExport(cl,
                            varlist = c("fun", "u1", "u2", "status1", "status2",
                                        "family", "x", "band", "kernel", 
                                        "optim_fun", "ieta"),
                            envir = environment())
    cveta <- parallel::parSapply(cl, X = seq_along(x), FUN = fun)
  }
  
  xind <- 1:length(x)
  nx <- length(x)
  obj <- SurvCopLocFun(u1 = u1[xind], u2 = u2[xind], 
                       status1 = status1[xind], status2 = status2[xind], 
                       family = family, x = cveta[xind], x0 = 0, eta = c(0,1),
                       wgt = rep(1, nx), degree = 1)
  cvll <- -obj$fn(c(0,1))
  if(!cveta_out) {
    return(cvll)
  } else {
    return(list(x = x, eta = cveta, loglik = cvll))
  }
}
