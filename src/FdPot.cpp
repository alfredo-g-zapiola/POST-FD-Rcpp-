#include "FdPot.h"
#include <assert.h>     /* assert */
#include <algorithm> // std::min_element
#include <cmath>

#include <Rcpp.h>
#include <omp.h>




using namespace fdpot;
using Rcpp::_;  // named placeholder, used to create an Rcpp List

template<>
ADdouble ORCT::cdf<ADdouble>(const ADdouble x) const{
   return 1 / (1 + CppAD::exp(- x * this->gamma));
};

Rcpp::List FdPot::fit(const arma::vec& y, const arma::mat & X_coeff,
                      const unsigned n_sols_){
  this->n_sols = n_sols_;
  // provide the tree class with the number of samples
#ifndef MYNDEBUG
  std::cout << "Received dataset with " << this->n_samples << "samples and " << 
    orct_ptr->n_labels << " different labels" << std::endl;
#else
  Rcpp::Rcout << "Received dataset with " << this->n_samples << " and " << 
    orct_ptr->n_labels << " different labels" << std::endl;
#endif 
  // what we are implicitly doing:
  // cast the coefficients Rcpp::NumericMatrix to an arma::mat (for vec multiplication)
  // avoid copies! https://stackoverflow.com/questions/57516726/converting-r-matrix-to-armamat-with-list
  // arma::mat X_cast = splines2::rmat2arma(X_coeff);

  // auto & [n1, n2] =  cannot do this in c++14
  // single class per leaf constraints
  this->n_constrs.first = orct_ptr->n_leaf_nodes;
  // at least one terminal node per class constraints
  this->n_constrs.second  = orct_ptr->n_labels;
#ifndef MYNDEBUG
  std::cout << "Number of single class per leaf constraints " << n_constrs.first << std::endl;
  std::cout << "Number of termianl node per class " << n_constrs.second << std::endl;
    std::cout << "Computing features" << std::endl;
#endif
#ifdef DEV 
  Rcpp::Rcout << "Computing features" << std::endl;
#endif
  // copy elision
  this->features = arma::mat(evalFd.compute_features(X_coeff, orct_ptr->n_feats));
  // this->features = X_coeff.t();
  // scale features
  this->scale_features(this->features);
#ifndef MYNDEBUG 
  std::cout << "Printing standardised features" << std::endl;
  for (unsigned p = 0; p < this->features.n_cols; p++)
    std::cout << this->features.col(p) << std::endl;
  std::cout << "Computing dissimilarity matrix" << std::endl;
#endif
#ifdef DEV 
  Rcpp::Rcout << "Computing dissimilarity matrix" << std::endl;
#endif
  this->dissim_matrix = this->evalFd.compute_dissim_matrix(X_coeff);
  
#ifndef MYNDEBUG 
  std::cout << "Setting up optimiser" << std::endl;
#endif
#ifdef DEV 
  Rcpp::Rcout << "Setting up optimiser" << std::endl;
#endif
  
  this->setup_optimiser(y, X_coeff);
    
#ifndef MYNDEBUG 
  std::cout << "Solving problems with different init points" << std::endl;
#endif
#ifdef DEV 
  Rcpp::Rcout << "Solving problems with different init points" << std::endl;
#endif
  
  auto res = this->solve_trees();
  
  return res;
}


OptimTraits::OptimFuns FdPot::create_mathematical_model(const arma::vec& y,
                                                        const arma::mat & X_coeff){
  
  
  // initialise what we retur  (recall std::function is a pointer wrapper)
  OptimTraits::OptimFuns f = nullptr;
  
  this->cost_func = [this, &y] (const VariantVarsT& vars) -> ADdouble {
    // I first create the visitor:
    
    
    ADdouble e_cost = 0.;
    // note this loop cannot be parallelised since it is the function that Ipopt will use
    for (unsigned i = 0; i < this->n_samples; i++){ // all samples
      for (unsigned leaf = orct_ptr->n_int_nodes; // first leaf comes after int_nodes
              leaf < orct_ptr->n_nodes; leaf++){  // all leafs
        
        ADdouble leaf_c = 0.;  // cost of the current leaf
        
        // first_class var in this leaf
        unsigned c_kt = orct_ptr->var_map(leaf);
        
        for (unsigned k = 0; k < orct_ptr->n_labels; k ++){
          // Sadly, I had to redefine the overloaded visitor at each iteration:
          auto visitor_add_leaf_c = OverloadedVisitor{
            [this, &leaf_c,
             &leaf, &y, &i, &k, &c_kt](const arma::vec& vars){
              leaf_c += orct_ptr->proba_fall_leaf<arma::vec>(
                features.row(i), vars, leaf) * this->miss_class_cost(y(i), k) * vars[c_kt++];
              // decided to use the ARMA_NO_DEBUG macro
            },
            [this, &leaf_c, &leaf, &y, &i, &k, &c_kt](
                const OptimTraits::ADvector& vars){
              leaf_c += orct_ptr->proba_fall_leaf<OptimTraits::ADvector>(
                this->features.row(i), vars, leaf) *\
                  this->miss_class_cost(y(i), k) * vars[c_kt++];
            }
          };
          
          
          //leaf_c += orct_ptr->proba_fall_leaf(this->features.row(i), vars, leaf) * 
          //  this->miss_class_cost(y(i), k) * vars[c_kt++];
          std::visit(visitor_add_leaf_c, vars);
          
        }
        e_cost += leaf_c;
      }
      
    }

    e_cost /= this->n_samples;
    return e_cost;
  };
  
  this->penalty_func = [this, &y] (const OptimTraits::ADvector& vars) -> ADdouble{
    // expected dissimilarity per leaf
    ADdouble e_diss = 0.;
    // the idx of the first leaf nodes is the number of interior nodes
    for (unsigned tau = orct_ptr->n_int_nodes; tau < orct_ptr->n_nodes; tau++){
      ADdouble leaf_d = 0.;  // leaf dissimilarity
      for (unsigned i = 0; i < this->n_samples; i++){
        for (unsigned j = i+1; j < this->n_samples; j++){
          leaf_d +=  orct_ptr->proba_fall_leaf(this->features.row(i), vars, tau) *
            orct_ptr->proba_fall_leaf(this->features.row(j), vars, tau) * 
            this->dissim_matrix(i, j);
        }
      }
      e_diss += leaf_d;
    }
    e_diss /= orct_ptr->n_leaf_nodes;
    return e_diss;
  };

  // obj_fn era qui    
  
  auto constr_single_class_pred = [this] (const OptimTraits::ADvector& vars,
                                          const unsigned leaf) -> ADdouble{
    ADdouble sum = 0.;
    unsigned first_idx = orct_ptr->var_map(leaf);
    for (unsigned k = 0; k < orct_ptr->n_labels; k++)
      sum += vars.at(first_idx++);  
    return sum;
  };
  
  auto constr_one_leaf_per_class = [this] (const OptimTraits::ADvector& vars,
                                     const unsigned k) -> ADdouble{
      ADdouble sum = 0.;
      for (unsigned leaf = orct_ptr->n_int_nodes; leaf < orct_ptr->n_nodes; leaf++)
        sum += vars.at(orct_ptr->var_map(leaf) + k);
      return sum; 
  };
  // prepare optim funs (functions passed by copy since else we get seg fault)
  f = [this, constr_single_class_pred, constr_one_leaf_per_class] 
      (OptimTraits::ADvector& fg, const OptimTraits::ADvector& x) -> void
    {
    // objective function
    // fg[0] = this->cost_func(x);
     fg[0] = this->obj_function(x);
     
  
    // single class prediction per leaf constraints
    unsigned i = 1;  // constraint index
    for (unsigned leaf = orct_ptr->n_int_nodes; leaf < orct_ptr->n_nodes; leaf++)
      fg[i++] = constr_single_class_pred(x, leaf);
    
    // at least one leaf per class constraints
    for (unsigned k = 0; k < orct_ptr->n_labels; k++)
      fg[i++] = constr_one_leaf_per_class(x, k);
    };
  
  return f;
  };

void FdPot::setup_optimiser(const arma::vec& y, const arma::mat & X_coeff){

  this->optimiser->create_variables(orct_ptr->n_vars,
                                    orct_ptr->n_leaf_nodes + orct_ptr->n_labels);
  
#ifndef MYNDEBUG
  std::cout<< "Created variables: n_feats " << orct_ptr->n_feats << " n int nodes: " <<
    orct_ptr->n_int_nodes << " n_labels " << orct_ptr->n_labels << " n_leaf nodes " <<
      orct_ptr->n_leaf_nodes << " yields n_vars " <<orct_ptr-> n_vars << 
        " and n_constraints " << orct_ptr->n_leaf_nodes + orct_ptr->n_labels << std::endl;
#endif
  
  // set up variables' bounds
  Dvector vars_lb(optimiser->n_vars), vars_ub(optimiser->n_vars);
  
  for (unsigned v = 0; v < (orct_ptr->n_int_nodes * (orct_ptr->n_feats + 1));
  v++){
    vars_lb[v] = -1.;
    vars_ub[v] = 1.;
  }
  for (unsigned v = orct_ptr->n_int_nodes * (orct_ptr->n_feats + 1);
       v < this->optimiser->n_vars; v++){
    vars_lb[v] = 0.;
    vars_ub[v] = 1.;
  }
  // setup constraints bounds
  std::vector<double> g_lb(n_constrs.first + n_constrs.second),
    g_ub(n_constrs.first + n_constrs.second);
  
  for (unsigned c = 0; c < n_constrs.first; c++){
    g_lb[c] = 1.; 
    g_ub[c] = 1.;
  }
  for (unsigned c = n_constrs.first; c < g_ub.size(); c++){
    g_lb[c] = 1.;
    g_ub[c] = 1.0e19; // no upper bound
  }
  this->optimiser->set_math_program(this->create_mathematical_model(y, X_coeff),
                                    std::move(vars_lb), std::move(vars_ub),
                                    std::move(g_lb), std::move(g_ub)
                                    );
}


void FdPot::initialise_vars(std::uint32_t seed, Dvector & vars){
  std::mt19937 engine{seed};
  dirichlet_distribution<std::mt19937> 
    dirichlet_distrib( std::vector<double>(orct_ptr->n_labels, 1.) );
  
  // fill the vector of variables
  // dirichlet for the leaf node variables (must sum to 1)
  for (unsigned leaf = orct_ptr->n_int_nodes; 
       leaf < orct_ptr->n_nodes; leaf++){
    
    unsigned first_idx = orct_ptr->var_map(leaf);
    // sample
    std::vector<double> sample = dirichlet_distrib(engine);
#ifndef MYNDEBUG
  //  std::cout << "Printing dirichlet sample for leaf " << leaf << " of\
    //sample of length " << sample.size() << std::endl;
#endif
#ifdef DEV
   // Rcpp::Rcout << "Printing dirichlet sample for leaf " << leaf << " of\
  //  sample of length " << sample.size() << std::endl;
#endif
    // now fill the leaf variables
    unsigned i = 0;
    for (unsigned k= first_idx; k < first_idx + orct_ptr->n_labels; k++){
            vars.at(k) = sample.at(i++);
    #ifndef MYNDEBUG
    //  std::cout << vars.at(k) << "\t";
    #endif
    }
    #ifndef MYNDEBUG 
    //std::cout<< std::endl;
    #endif
  }
  // uniform for all other parameters
  std::uniform_real_distribution<> uniform_distrib(0, 1);  // set up distrib
  
  for (unsigned node = 0; node < orct_ptr->n_int_nodes; node++){
    unsigned first_idx = orct_ptr->var_map(node);
    
  #ifndef MYNDEBUG
   // std::cout << "Printing variables for node " << node << " whose\
    //starting index in the vars vector is " << first_idx << std::endl;
  #endif
    // we have n_feats + 1 variables per interior node
    for (unsigned j = first_idx; j < first_idx + (orct_ptr->n_feats + 1); j++){
      vars.at(j) = uniform_distrib(engine);
      #ifndef MYNDEBUG
       // std::cout<< vars.at(j) << "\t";
      #endif
    }
    #ifndef MYNDEBUG
     // std::cout << std::endl;
    #endif
  }
}
  
Rcpp::List FdPot::solve_trees(void){
  // create the random seeds
  std::vector<unsigned> seeds(this->n_sols);  // setup seeds vector
  for (unsigned i = 0; i < this->n_sols; i++)
    seeds.at(i) = this->seed + i;
  std::seed_seq seq({this->seed});
  seq.generate(seeds.begin(), seeds.end());
  
  // prepare results variable
  FdPotResults results(this->n_sols, orct_ptr->n_vars);

  std::vector<OptimHandler> optimhandlers(n_sols);
  
  #pragma omp parallel for
  for (unsigned m = 0; m < n_sols; m++){
    optimhandlers[m] = OptimHandler(*(this->optimiser));
    // initialise variables with current seed
   this->initialise_vars(seeds.at(m), optimhandlers[m].variables);
  }
  
#ifdef _OPENMP
  
  int threads = omp_get_num_threads();
  std::cout << "open mp active. n threads " << threads << std::endl;
  std::cout << "Runinng in parallel: " << std::endl;
#endif
  
  // #pragma omp parallel for shared(optimhandlers, results)
  for (unsigned m = 0; m < n_sols; m++){
    // initialise the current optimisation handler
    // NB the copy constructor is used since the structure is the same
    // for std::move another variable would have to be used; the gained efficiency
    // was not worth the "ugly" code in my opinion
    auto& cur_optim_hdler =  optimhandlers[m];  
    

    // perform the optimisation
    try{
      cur_optim_hdler.solve();
    }
    catch (const std::runtime_error&e){
      Rcpp::Rcout << e.what() << std::endl;
      continue;
    }
    
    results.all_variables.col(m) = std::move(cur_optim_hdler.solution.x);

    auto vars_ad = OptimTraits::ADvector(results.all_variables.col(m).cbegin(),
                                         results.all_variables.col(m).cend());
    
    //results.cost_func_vals(m) = CppAD::Value(this->cost_func(vars_ad));
    results.cost_func_vals(m) = CppAD::Value(this->cost_func(results.all_variables.col(m)));
      //CppAD::Value(*(this->cost_func)(cur_optim_hdler.variables));

    results.penalty_func_vals(m) = \
      CppAD::Value(this->penalty_func(vars_ad));
    
    results.obj_func_vals(m) =  cur_optim_hdler.solution.obj_value;
    Rcpp::checkUserInterrupt();  // check if the user has clicked stop

}  // end for
  //#ifndef MYNDEBUG
      //std::cout << "Solution " << m << " with value " <<
     //   cur_optim_hdler.solution.obj_value <<std::endl;
  //#endif
  unsigned best_idx = std::distance(results.obj_func_vals.cbegin(), 
                                    std::min_element(
                                       results.obj_func_vals.cbegin(),
                                       results.obj_func_vals.cend()
                                                    )
                                      );
#ifndef MYNDEBUG
  
  std::cout << "Sending out solutions list "  << std::endl;
  std::cout << "Depth: " <<  orct_ptr->depth << std::endl;
  std::cout << "n_feats: " <<  orct_ptr->n_feats << std::endl;
  std::cout << "n_labels: " <<  orct_ptr->n_labels << std::endl;
  std::cout << "best_idx: " <<  best_idx << std::endl;
#endif
    
  return Rcpp::List::create(
    _("depth") = orct_ptr->depth,
    _("n_feats") = orct_ptr->n_feats,
    _("n_labels") = orct_ptr->n_labels,
    _("best_tree_idx") =  best_idx,
    _("obj_func_vals") = results.obj_func_vals,
    _("cost_func_vals") = results.cost_func_vals,
    _("penalty_func_vals") = results.penalty_func_vals,
    _("all_variables") = results.all_variables,
    _("best_variables") = results.all_variables.col(best_idx)
    
    // _("obj_func_vals") = std::move(results.best_variables)
    //_["cost_func_vals"] = arma::vec(this->n_sols) ,
   // _["penalty_func_vals"] = arma::vec(this->n_sols),
  );
};

void FdPot::scale_features(arma::mat& features){
    #ifndef MYNDEBUG
      std::cout << "Scaling features" << std::endl;
    #endif
    #ifdef DEV
        Rcpp::Rcout << "Scaling features" << std::endl;
    #endif
    // perform MinMax scaling
    for (unsigned p = 0; p < features.n_cols; p++){
      double colmax = features.col(p).max();
      double colmin = features.col(p).min();
      
      features.col(p) = (features.col(p) - colmin) / (colmax - colmin);
    }
};



