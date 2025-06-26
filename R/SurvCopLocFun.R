#' Create a \pkg{TMB} local likelihood function for survival data.
#'
#' Wraps a call to [TMB::MakeADFun()].
#'
#' @template param-u1
#' @template param-u2
#' @param status1 Vector of censoring indicators for the first variable.
#' @param status2 Vector of censoring indicators for the second variable.
#' @template param-family
#' @template param-x
#' @param x0 Scalar covariate value at which to evaluate the local likelihood.  Does not have to be a subset of `x`.
#' @param wgt Vector of positive kernel weights.
#' @template param-degree
#' @param eta Value of the copula dependence parameter.  Scalar or vector of length two, depending on whether `degree` is 0 or 1.
#' @param rm.zero Logical indicating whether zero values in `u1` and `u2` should be removed when calculating the local likelihood. If `rm.zero = TRUE`, the rows with `u1=0` or `u2=0` are removed from the data, if `rm.zero = FALSE`, then a small quantity is added to `u1` and `u2` to avoid numerical issues. 
#' @return A list as returned by a call to [TMB::MakeADFun()].  In particular, this contains elements `fun` and `gr` for the *negative* local likelihood and its gradient with respect to `eta`.
#' @export
SurvCopLocFun <- function(u1, u2, status1, status2, 
                          family, x, x0, wgt, 
                          degree = 1, eta, 
                          rm.zero = TRUE) {
  if(!family %in% 3:5) {
    stop("Unsupported copula family (must be integer between 3-5).")
  }
  
  # subset data 
  wpos <- wgt > 0 # index of positive weights
  u1 <- u1[wpos]
  u2 <- u2[wpos]
  status1 <- status1[wpos]
  status2 <- status2[wpos]
  wgt <- wgt[wpos]
  x <- x[wpos]
  
  # handle zeros in u1 and u2
  iz <- which(u1==0 | u2==0)
  if(length(iz) > 0){
  if(rm.zero == TRUE){
    u1 <- u1[-iz]
    u2 <- u2[-iz] 
    status1 <- status1[-iz]
    status2 <- status2[-iz]
    x <- x[-iz]
    wgt <- wgt[-iz]
  }else{
    u1 <- ifelse(u1==0, u1 + runif(1, min = 0, max = 10^-8), u1)
    u2 <- ifelse(u2==0, u2 + runif(1, min = 0, max = 10^-8), u2)
  }} 
  

  # create TMB function
  # censoring groups 
  delta1 <- (1-status1)*(1-status2) # both censored
  delta2 <- status1*(1-status2) # first uncensored, second censored
  delta3 <- (1-status1)*status2 # first censored, second uncensored
  delta4 <- status1*status2 # both uncensored
  ix <- c(which(delta4==1), which(delta3==1), which(delta2==1), which(delta1==1))
  # data input
  data <- list(model = "SurvLikelihood",
               u1 = u1[ix], u2 = u2[ix],
               cen_start = c(0, sum(delta4), sum(delta4)+sum(delta3),
                             sum(delta4)+sum(delta3)+sum(delta2)), 
               cen_length = c(sum(delta4), sum(delta3), sum(delta2), sum(delta1)),
               wgt = wgt[ix], xc = x[ix]-x0,
               family = family)
  parameters <- list(beta = eta)
  # convert degree to TMB::map
  ## degree <- .format_degree(degree)
  if(!degree %in% 0:1) stop("degree must be 0 or 1.")
  map <- list(beta = factor(c(1, 2)))
  if(degree == 0) {
    parameters$beta[2] <- 0
    map$beta[2] <- NA
  }
  TMB::MakeADFun(
    data = data,
    parameters = parameters,
    map = map,
    DLL = "LocalCop_TMBExports",
    silent = TRUE
  )
}


