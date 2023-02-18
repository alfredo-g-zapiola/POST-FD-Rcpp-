// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <iostream>




// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <splines2Armadillo.h>
#include "FdPot.h"
#include "helpers.h"
  
#include <functional>
#include <dlfcn.h>
  

using Rcpp::_; //aka Named (used to create an Rcpp::List)
using namespace fdpot;

//' Build and fit an FD-classification penalised tree
//' 
//' @description instantiates and fits a Functional Data Penalised Optimial Randomised Decision Tree
//' @note The s3 class p.fdorct provides an interface for this method 
//' @param y labels vector. Note they have to be integers starting from 0
//' @param X_coeffs p x n matrix with the coefficients fitted in the smoothing (see the examples of the library), where p is the number of coeffients for func. datum and n the sample size
//' @param X_argvals vector of lower and upper bounds of the domain
//' @param X_basis_df the degrees of freedom of the basis
//' @param X_basis_degree the degree of the (b-spline) basis
//' @param basis_type the basis type string, BSpline is the only one currently supported
//' @param depth the tree depth. Make 
//' @param alpha the hyperparameter for the penalty in the objective function. Higher alpha, higher weight for the penalty
//' @þaram similarity method: the method to compute similarity between two functions (by default, L2 norm)
//' @param n_feats how features to use at each node of the tree to perform a split. An equal-length partition of the size of n_feats is created; each feat is the integral of the func. datum in a set that is part of the pariition
//' @param n_solve how many different trees to fit starting from different init points
//' @param gamma the randomisation factor for the ORCT, best kept default
//' @param seed random seed for reproducibility
// [[Rcpp::export]]
Rcpp::List pFdorct_Rcpp(const arma::vec & y, 
                        const arma::mat&  X_coeffs,
                        const Rcpp::NumericVector & X_argvals,
                        int X_basis_df,
                        int X_basis_degree,
                        const Rcpp::String & basis_type = "BSpline",
                        int depth = 2,
                        double alpha = .1,
                        Rcpp::String similarity_method = "d0.L2",
                        unsigned n_feats = 10,
                        int n_solve = 20 ,
                        double gamma = 512.,
                        long int seed = 41703192
){
   //1 Basis object
   // Call template class with basis, params

#ifdef DEV
   Rcpp::Rcout << "DEV defined" << std::endl;
#endif
 // Call template class with basis, params
// #define MYNDEBUG
#ifdef MYNDEBUG
   Rcpp::Rcout << "MYNDEBUG defined" << std::endl;
#endif
#ifndef MYNDEBUG
   Rcpp::Rcout << "MYNDEBUG not defined" << std::endl;
#endif
   
#ifndef MYNDEBUG 
  std::cout << "Obtaining basis object" << std::endl;
  std::cout << "Boundaries: " << X_argvals[0] << " and "\
    << X_argvals[X_argvals.size()-1] << std::endl;
#endif
  // obtain the boundary knots
#ifdef DEV
  Rcpp::Rcout << "Obtaining boundary knots" << std::endl;
#endif
  arma::vec boundary_knots{ X_argvals[0], X_argvals[X_argvals.size()-1] };
  
#ifndef MYNDEBUG
  std::cout << "Printing boundary knots: \n" << boundary_knots << std::endl;
#endif
#ifdef DEV
  Rcpp::Rcout << "Constructing basis" << std::endl;
#endif
  // now construct the basis
  auto basis = splines2::BSpline(X_argvals, X_basis_df, X_basis_degree,
                                 boundary_knots);
  // Report/PresentationRaccontare dynamic loading (storia)
  
  // Validate dimensions
  if (y.size() != X_coeffs.n_cols){
  #ifndef MYNDEBUG
    std::cerr << "Error when working without MYNDEBUG" << std::endl;
    std::cerr<<"Number of samples (cols in the coefficients matrix) must conform\
                 the number of labels: provided " << y.size() << " labels"      \
                                                  << " and " << X_coeffs.n_cols << " samples" << std::endl;
               //   std::exit(1);
                 
                 Rcpp::stop("Number of rows in the coefficients matrix must\
                 conform to the number of labels");
  #else
                 Rcpp::stop("Number of rows in the coefficients matrix must\
                 conform to the number of labels");
  #endif 
  }

  unsigned n_samples = X_coeffs.n_cols;
  unsigned n_labels = arma::vec(arma::unique(y)).n_rows ; 

  // Validate depths:
  if (not (depth > 0)){
    Rcpp::stop("depth must be at least 1");
  }
  // recall cum_tree_sum counts the number of nodes WITHOUT including the current one
  unsigned n_leaf_nodes = helpers::cum_tree_sum(depth+1) - \
                          helpers::cum_tree_sum(depth);
#ifndef MYNDEBUG
  std::cout << "number of leaf nodes: " << n_leaf_nodes << std::endl;
  std::cout << "number of labels: " << n_labels << std::endl;
#endif 
#ifdef DEV
  Rcpp::Rcout <<"number of leaf nodes: " << n_leaf_nodes << std::endl;
  Rcpp::Rcout <<  "number of labels: " << n_labels << std::endl;
#endif
  if(n_leaf_nodes < n_labels){
    #ifndef MYNDEBUG
    std::cerr << "Number of leaf nodes must be >= the number of labels,\
               increase the depth" << std::endl;
    std::exit(1);
    #else
    Rcpp::stop("Number of leaf nodes must be >= the number of labels,\
               increase the depth" );
    #endif
  }
  
      

 
 // use the Communicator (like in the design pattern)
 // FunFact I became fond of such DP in the last challenge for the PACS course   
  #ifdef DEV
  Rcpp::Rcout << "Instantiating tree" << std::endl;
  #endif 
  FdPot tree = FdPot(std::move(basis), n_labels, n_samples, n_feats, depth, alpha,
                    seed, gamma);
  #ifdef DEV
  Rcpp::Rcout << "Fitting tree" << std::endl;
  #endif
  
  Rcpp::List fit_results = tree.fit(y, X_coeffs, n_solve);
  Rcpp::List fitted_tree_specs = Rcpp::List::create(
    _("fit_results") = fit_results,
    _("n_samples") = n_samples,
    _("alpha") = alpha,
    _("X_argvals") = X_argvals,
    _("X_basis_df") = X_basis_df,
    _("X_basis_degree") = X_basis_degree,
    _("boundary_knots") = boundary_knots,
    _("gamma") = gamma,
    _("seed") = seed
  );
    
  return fitted_tree_specs;
  };

// [[Rcpp::export]]
Rcpp::List predict_FdPot_Rcpp(const Rcpp::List& fitted_tree,
                        const arma::mat& X_coefs){
  Rcpp::List fit_results = Rcpp::as<Rcpp::List>(fitted_tree["fit_results"]);
  // to obtain the features, need the basis
  auto basis = splines2::BSpline(
    Rcpp::as<Rcpp::NumericVector>(fitted_tree["X_argvals"]), 
    Rcpp::as<int>(fitted_tree["X_basis_df"]), 
    Rcpp::as<int>(fitted_tree["X_basis_degree"]),
    Rcpp::as<arma::vec>(fitted_tree["boundary_knots"])
  );
  
  FdHandler<BasisEnum::BSPLINE> fd_handler(std::move(basis));
  auto feats = arma::mat( fd_handler.compute_features(X_coefs, 
                                                      fit_results("n_labels"))
                          );
  auto vars = Rcpp::as<arma::vec>(fit_results["best_variables"]);// TODO
  
  // scale features
  FdPot::scale_features(feats);
#ifndef MYNDEBUG
  std::cout << "Creating ORCT with following params" << std::endl;
  std::cout << "Depth " << Rcpp::as<int>(fit_results["depth"]) << std::endl;
  std::cout << "N_feats " <<   Rcpp::as<int>(fit_results["n_feats"]) << std::endl;
  std::cout << "n labels: " <<  Rcpp::as<int>(fit_results["n_labels"]) << std::endl;
#endif 
  ORCT tree = ORCT(fit_results["depth"], fit_results["n_feats"], 
                   fit_results["n_labels"], fitted_tree["gamma"]
                     );

  return tree.predict(feats, vars);
}
