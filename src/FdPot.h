#ifndef p_fdorct
#define p_fdorct
 // class enum similarity method todo
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
    
    inline double miss_class_cost(const unsigned label, const unsigned predicted) const {
      return (label == predicted) ? 0. : this->missclaf_cost;
    }
    
    using VariantVarsT = std::variant<OptimTraits::ADvector, arma::vec>;
    
    arma::mat dissim_matrix;
    arma::mat features;
  private:
    FdHandler<BasisEnum::BSPLINE> evalFd;  // bridge to value functional data.
    std::unique_ptr<ORCT> orct_ptr = nullptr;
    std::unique_ptr<OptimHandler> optimiser = std::make_unique<OptimHandler>();
  	
  	std::pair<unsigned, unsigned> n_constrs = std::make_pair(0,0); // updated in the fit method
  	unsigned long seed;
  	double alpha;
  	unsigned n_sols = 0u;
  	double missclaf_cost{0.5}; // misclassification cost, this number was used in the experiments by Blaquero et al.
  	unsigned n_samples = 0u;
  	// 1Rcpp::String similarity_method; // TODO.
  	//////////////////////////////////////////////////////////

  
  	std::function<ADdouble (const VariantVarsT&)> cost_func = nullptr;
  	std::function<ADdouble(const OptimTraits::ADvector&)> penalty_func = nullptr;
    //std::function<ADdouble(const std::variant<OptimTraits::ADvector, 
      //                    arma::vec> & vars)> a = nullptr; //  obj_function = [this](const OptimTraits::ADvector& vars) -> ADdouble{
  	 // return this->cost_func(vars) + this->alpha * this->penalty_func(vars);
  // 	};
  	

  	
  	Rcpp::List solve_trees(void);
  

  	
  	OptimTraits::OptimFuns create_mathematical_model(const arma::vec& y, const arma::mat & X_coeff); 
  	
  	void setup_optimiser(const arma::vec& y, const arma::mat & X_coeff);
  	
  	void initialise_vars(std::uint32_t seed, OptimHandler::Dvector & vars);
  	
  	
  	
  }; // class FdPot
    
// Here I pasted Prof. Formaggia's variadic template for overloading a Visitor

// This part is only for the nerds of you (skip if you are not interested)
// This is a class with which you can overload an arbitrary number of functors
template <class... Ts> struct OverloadedVisitor : Ts...
{
  using Ts::operator()...;
};
// explicit deduction guide (not explained at lecture: it allows to create an
// overloaded object just using the constructor. Very nice, but it was too much
// for the course....). Moreover it is not necessary since C+=20.
template <class... Ts> OverloadedVisitor(Ts...) -> OverloadedVisitor<Ts...>; 
// End nerdish part
    
//  The other options was a custom visitor: see https://www.cppstories.com/2018/09/visit-variants/
// This would have required a constructor and constructing at each iteration so not worth it
//struct Visitor_add_leaf_cost{
//   
//   
//   void operator()(const arma::vec& vars, ADdouble& leaf_c,
//                 const std::unique_ptr<ORCT>& orct_ptr,
//                 const FdPot* & that ,
//                 const unsigned leaf,
//                 const arma::vec& y,
//                 const unsigned i,
//                 const unsigned k,
//                 const unsigned c_kt){
//     leaf_c += orct_ptr->proba_fall_leaf<arma::vec>(that->features.row(i), vars, leaf) * 
//       that->miss_class_cost(y(i), k) * vars[c_kt];
//   }
//   
//   void operator()(const OptimTraits::ADvector& vars, ADdouble& leaf_c, 
//                 const std::unique_ptr<ORCT>& orct_ptr,
//                 const FdPot* &that,
//                 const unsigned leaf,
//                 const arma::vec& y,
//                 const unsigned i,
//                 const unsigned k,
//                 const unsigned c_kt){
//     leaf_c += orct_ptr->proba_fall_leaf<OptimTraits::ADvector>(that->features.row(i), vars, leaf) * 
//       that->miss_class_cost(y(i), k) * vars[c_kt];}
// }

 }; // namespace fdpot

#endif // p_fdorct
