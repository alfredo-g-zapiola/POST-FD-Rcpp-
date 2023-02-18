#ifndef p_fdorct
#define p_fdorct
 // TODO size_t
 // class enum similarity method
 // TODO documentation (doxygen)s
#include <cmath>
#include <tuple>
#include <random>
#include <variant>
#include <omp.h>
 
#include "BasisObj.h"
#include "OptimHandler.h"
#include "ORCT.h"
#include "helpers.h"
#include "FdPotSupport.h"
#include "dirichlet.h"
 
 namespace fdpot{
 // add the specialisation for ADdouble (we may need ORCT in another context)

struct VisitorProbaGoLeft{
  inline ADdouble operator()(const OptimTraits::ADvector& vars, 
                           const std::unique_ptr<ORCT> orct_ptr&)const {
      
      return orct_ptr->proba_fall_leaf(feats, vars, tau)
    
    
  }
};

class FdPot: public OptimTraits{
public:
    FdPot(splines2::BSpline&& basis_,
          const unsigned n_labels,
          const unsigned n_samples_,
          const unsigned n_feats,
          const unsigned int depth_,
          const double alpha_,
          const unsigned long int seed_ = 22200337,
          const double gamma_=512) : 
    orct_ptr{std::make_unique<ORCT>(depth_, n_feats, n_labels, gamma_)},
    evalFd{std::move(basis_)},
    n_samples(n_samples_),
    alpha{alpha_}, 
    seed{seed_}
    {};
    
    Rcpp::List fit(const arma::vec& y, const arma::mat& X_coeff, const unsigned n_sols = 20);
    
    static void scale_features(arma::mat&);
    
    using VariantVarsT = std::variant<OptimTraits::ADvector, arma::vec>;
    
  
  private:
    FdHandler<BasisEnum::BSPLINE> evalFd;  // bridge to value functional data.
    std::unique_ptr<ORCT> orct_ptr = nullptr;
    std::unique_ptr<OptimHandler> optimiser = std::make_unique<OptimHandler>();
  	
  	std::pair<unsigned, unsigned> n_constrs = std::make_pair(0,0); // TODO
  	unsigned long seed = 100ul;
  	double alpha;
  	unsigned n_sols = 0u;
  	double missclaf_cost{0.5}; // misclassification cost, this number was used in the experiments by Blaquero et al.
  	unsigned n_samples = 0u;
  	// 1Rcpp::String similarity_method; // TODO.
  	//////////////////////////////////////////////////////////
  	arma::mat dissim_matrix;
  	arma::mat features;
  
  	std::function<ADdouble (const VariantVarsT &)> cost_func = nullptr;
  	std::function<ADdouble(const OptimTraits::ADvector&)> penalty_func = nullptr;
    //std::function<ADdouble(const std::variant<OptimTraits::ADvector, 
      //                    arma::vec> & vars)> a = nullptr; //  obj_function = [this](const OptimTraits::ADvector& vars) -> ADdouble{
  	 // return this->cost_func(vars) + this->alpha * this->penalty_func(vars);
  // 	};
  	
  	
  		///////////////////////////////////////////
  	
  	Rcpp::List solve_trees(void);
  
  	inline double miss_class_cost(const unsigned label, const unsigned predicted){
  	  return (label == predicted) ? 0. : this->missclaf_cost;
  	}
  	
  	OptimTraits::OptimFuns create_mathematical_model(const arma::vec& y, const arma::mat & X_coeff); 
  	
  	void setup_optimiser(const arma::vec& y, const arma::mat & X_coeff);
  	
  	void initialise_vars(std::uint32_t seed, OptimHandler::Dvector & vars);
  	
  	
  	
  }; // class FdPot
    
  
 }; // namespace fdpot

#endif // p_fdorct
