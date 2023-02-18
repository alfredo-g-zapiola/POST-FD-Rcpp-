#ifndef FDPOT_SUPPORT_HEADER
#define FDPOT_SUPPORT_HEADER
#include "BasisObj.h"

struct FdPotOptions{
  splines2::BSpline basis_;
  
  
};

struct FdPotResults{
  FdPotResults(const unsigned n_sols, const unsigned n_vars):
    obj_func_vals(arma::vec(n_sols)), cost_func_vals(arma::vec(n_sols)),
    penalty_func_vals(arma::vec(n_sols)), 
    all_variables(arma::mat(n_vars, n_sols))
  {};
  
  arma::vec obj_func_vals, cost_func_vals, penalty_func_vals, best_variables;
  arma::mat  all_variables;
  
  
};

#endif 