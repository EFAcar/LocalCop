// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_LocalCop_TMBExports
#include <TMB.hpp>
#include "ClaytonCDF.hpp"
#include "ClaytonNLL.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "ClaytonCDF") {
    return ClaytonCDF(this);
  } else if(model == "ClaytonNLL") {
    return ClaytonNLL(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
