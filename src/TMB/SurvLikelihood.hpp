/// @file SurvLikelihood.hpp

#include "LocalCop/frank.hpp"
#include "LocalCop/gumbel.hpp"
#include "LocalCop/clayton.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// Weighted conditional copula loglikelihood for survival data.
///
/// @param[in] u1 First response vector of length `n_obs`.
/// @param[in] u2 Second response vector of length `n_obs`.
/// @param[in] cen_start Vector of length 4 giving the start locations for the uncensored, first censored, second censored, and both censored observations.
/// @param[in] cen_length Vector of length 4 giving the length of the uncensored, first censored, second censored, and both censored observations.
/// @param[in] wgt Weight vector of length `n_obs`.
/// @param[in] xc Centered covariate vector of length `n_obs`.
/// @param[in] family Integer specifying the copula family (Clayton = 3, Gumbel = 4, Frank = 5).
/// @param[in] beta Parameter vector of length 2.
template<class Type>
Type SurvLikelihood(objective_function<Type> *obj) {
  DATA_VECTOR(u1); // first response vector
  DATA_VECTOR(u2); // second response vector
  DATA_IVECTOR(cen_start); // a vector of length 4 giving the start locations for uncensored, first censored, second censored, and both censored. 
  DATA_IVECTOR(cen_length); // a vector of length 4 giving the size of the uncensored, first censored, second censored, and both censored. 
  DATA_VECTOR(wgt); // weights
  DATA_VECTOR(xc); // centered covariates, i.e., X - x
  DATA_INTEGER(family); // copula family: 3-5.
  PARAMETER_VECTOR(beta); // dependence parameter: eta = beta[0] + beta[1] * xc
  
  Type nll = 0.0;
  vector<Type> eta = beta(0) + beta(1) * xc;
  Type res = Type(0.0);
  if(family == 3) {
    // Clayton copula
    vector<Type> theta = eta.exp();
    // uncensored
    int j=0;
    if(cen_length(j) > 0) {
      vector<Type> _wgt = wgt.segment(cen_start(j), cen_length(j));
      vector<Type> _u1 = u1.segment(cen_start(j), cen_length(j));
      vector<Type> _u2 = u2.segment(cen_start(j), cen_length(j));
      vector<Type> _theta = theta.segment(cen_start(j), cen_length(j));
      res += (_wgt * LocalCop::dclayton(_u1, _u2, _theta, 1)).sum();
      // res += (wgt.segment(cen_start(j), cen_length(j)) * 
      //   LocalCop::dclayton(u1.segment(cen_start(j), cen_length(j)), 
      //                      u2.segment(cen_start(j), cen_length(j)), 
      //                      theta.segment(cen_start(j), cen_length(j)), 1)).sum();
      }
    // first censored
    j=1;
    if(cen_length(j) > 0) {
      vector<Type> _wgt = wgt.segment(cen_start(j), cen_length(j));
      vector<Type> _u1 = u1.segment(cen_start(j), cen_length(j));
      vector<Type> _u2 = u2.segment(cen_start(j), cen_length(j));
      vector<Type> _theta = theta.segment(cen_start(j), cen_length(j));
      res += (_wgt * LocalCop::hclayton(_u1, _u2, _theta, 1)).sum();
      // res += (wgt.segment(cen_start(j), cen_length(j)) *
      //   LocalCop::hclayton(u1.segment(cen_start(j), cen_length(j)),
      //                      u2.segment(cen_start(j), cen_length(j)),
      //                      theta.segment(cen_start(j), cen_length(j)), 1)).sum();
      }
    // second censored
    j=2;
    if(cen_length(j) > 0) {
      vector<Type> _wgt = wgt.segment(cen_start(j), cen_length(j));
      vector<Type> _u1 = u1.segment(cen_start(j), cen_length(j));
      vector<Type> _u2 = u2.segment(cen_start(j), cen_length(j));
      vector<Type> _theta = theta.segment(cen_start(j), cen_length(j));
      res += (_wgt * LocalCop::hclayton(_u2, _u1, _theta, 1)).sum();
      // res += (wgt.segment(cen_start(j), cen_length(j)) * 
      //   LocalCop::hclayton(u2.segment(cen_start(j), cen_length(j)), 
      //                      u1.segment(cen_start(j), cen_length(j)), 
      //                      theta.segment(cen_start(j), cen_length(j)), 1)).sum();
    }
    // both censored
    j=3;
    if(cen_length(j) > 0) {
      vector<Type> _wgt = wgt.segment(cen_start(j), cen_length(j));
      vector<Type> _u1 = u1.segment(cen_start(j), cen_length(j));
      vector<Type> _u2 = u2.segment(cen_start(j), cen_length(j));
      vector<Type> _theta = theta.segment(cen_start(j), cen_length(j));
      res += (_wgt * LocalCop::pclayton(_u1, _u2, _theta, 1)).sum();
      // res += (wgt.segment(cen_start(j), cen_length(j)) * 
      //   LocalCop::pclayton(u1.segment(cen_start(j), cen_length(j)), 
      //                      u2.segment(cen_start(j), cen_length(j)), 
      //                      theta.segment(cen_start(j), cen_length(j)), 1)).sum();
      }
    // } else if(family == 4) {
    //   // Gumbel copula
    // vector<Type> theta = 1.0 + eta.exp();
    //   // uncensored
    //   int j=0;
    //   if(cen_length(j) > 0) {
    //     res += wgt.segment(cen_start(j), cen_length(j)) * LocalCop::dgumbel(
    //       u1.segment(cen_start(j), cen_length(j)), 
    //       u2.segment(cen_start(j), cen_length(j)), 
    //       theta.segment(cen_start(j), cen_length(j)), 1).sum();
    //   }
    //   // first censored
    //   j=1;
    //   if(cen_length(j) > 0) {
    //     res += wgt.segment(cen_start(j), cen_length(j)) * LocalCop::hgumbel(
    //       u1.segment(cen_start(j), cen_length(j)), 
    //       u2.segment(cen_start(j), cen_length(j)), 
    //       theta.segment(cen_start(j), cen_length(j)), 1).sum();
    //   }
    //   // second censored
    //   j=2;
    //   if(cen_length(j) > 0) {
    //     res += wgt.segment(cen_start(j), cen_length(j)) * LocalCop::hgumbel(
    //       u2.segment(cen_start(j), cen_length(j)), 
    //       u1.segment(cen_start(j), cen_length(j)), 
    //       theta.segment(cen_start(j), cen_length(j)), 1).sum();
    //   }
    //   // both censored
    //   int j=3;
    //   if(cen_length(j) > 0) {
    //     res += wgt.segment(cen_start(j), cen_length(j)) * LocalCop::pgumbel(
    //       u1.segment(cen_start(j), cen_length(j)), 
    //       u2.segment(cen_start(j), cen_length(j)), 
    //       theta.segment(cen_start(j), cen_length(j)), 1).sum();
    //   }
    // } else if(family == 5) {
    //   // Frank copula
    //   vector<Type> theta = eta;
    //   // uncensored
    //   int j=0;
    //   if(cen_length(j) > 0) {
    //     res += wgt.segment(cen_start(j), cen_length(j)) * LocalCop::dfrank(
    //       u1.segment(cen_start(j), cen_length(j)), 
    //       u2.segment(cen_start(j), cen_length(j)), 
    //       theta.segment(cen_start(j), cen_length(j)), 1).sum();
    //   }
    //   // first censored
    //   int j=1;
    //   if(cen_length(j) > 0) {
    //     res += wgt.segment(cen_start(j), cen_length(j)) * LocalCop::hfrank(
    //       u1.segment(cen_start(j), cen_length(j)), 
    //       u2.segment(cen_start(j), cen_length(j)), 
    //       theta.segment(cen_start(j), cen_length(j)), 1).sum();
    //   }
    //   // second censored
    //   int j=2;
    //   if(cen_length(j) > 0) {
    //     res += wgt.segment(cen_start(j), cen_length(j)) * LocalCop::hfrank(
    //       u2.segment(cen_start(j), cen_length(j)), 
    //       u1.segment(cen_start(j), cen_length(j)), 
    //       theta.segment(cen_start(j), cen_length(j)), 1).sum();
    //     }
    //   }
    //   // both censored
    //   int j=3;
    //   if(cen_length(j) > 0) {
    //     res += wgt.segment(cen_start(j), cen_length(j)) * LocalCop::pfrank(
    //       u1.segment(cen_start(j), cen_length(j)), 
    //       u2.segment(cen_start(j), cen_length(j)), 
    //       theta.segment(cen_start(j), cen_length(j)), 1).sum();
    //   }
    } else {
      Rf_error("Unknown copula family.");
    }
  nll = -res;
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
