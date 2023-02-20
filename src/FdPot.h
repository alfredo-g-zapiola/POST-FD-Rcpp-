#ifndef p_fdorct
#define p_fdorct

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
/*! @brief Constructor
	Initialises the pointer to the ORCT; the FdHandler
@param n_labels the number of different classes in the classification problem.
@param n_samples_ the number of statistical units
@param n_feats the number of features to calculate for each functional datum
@param depth_ the tree depth
@param alpha_ the penalisation weight
@param seed_ the seed for the different initialisation points
@param gamma_ randomisation factor, best left unchanged

*/
    FdPot(splines2::BSpline&& basis_,
          const unsigned n_labels,
          const unsigned n_samples_,
          const unsigned n_feats,
          const unsigned int depth_,
          const double alpha_,
          const unsigned long int seed_ = 22200337,
          const double gamma_=512.) : 
    orct_ptr{std::make_unique<ORCT>(depth_, n_feats, n_labels, gamma_)},
    evalFd{std::move(basis_)},
    n_samples(n_samples_),
    alpha{alpha_}, 
    seed{seed_}
    {};
    
    /*! @brief Calls different methods to orchestrate fitting
    @param y the labels vector
    @param X_coeff the coefficients matrix of the smoothing
    @param n_sols the number of solution
    @return an Rcpp::List with the fitting results
    */
    Rcpp::List fit(const arma::vec& y, const arma::mat& X_coeff, const unsigned n_sols = 20);
    
    /*! @brief Scales the features between 0 and 1 
    It is a MinMax Scaler, moving everything to [0,1]
    @param feats  the features matrix
    Note it transforms them in-place (it is void)
    */
    static void scale_features(arma::mat& feats);
     /*! @brief Calculates misclassificaton cost
     @param label the actual albel
     @param predicted the predicted label
    */   
    inline double miss_class_cost(const unsigned label, const unsigned predicted) const {
      return (label == predicted) ? 0. : this->missclaf_cost;
    }
     /*! @brief Variant type for the variables
     
     It can store either an ADvector, as used by the OptimHandler, or 
     and arma::vec, when passing the solutions of the optimisation
    */   
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
  	// 1Rcpp::String similarity_method; // for the future
  	//////////////////////////////////////////////////////////

      /*! @brief Pointer to the cost function
      @note it is updated runtime
    	@param vars the variables
    	*/
  	std::function<ADdouble (const VariantVarsT& vars)> cost_func = nullptr;
      /*! @brief the penalty function pointer
      Updated runtime
    */
  	std::function<ADdouble(const OptimTraits::ADvector&)> penalty_func = nullptr;
      /*! @brief Objective function
      The linear combination of the expected misclassification cost and the 
      dissimilarity penalisation.
    */
  	std::function<ADdouble(const OptimTraits::ADvector&)>  obj_function = \
  	  [this](const OptimTraits::ADvector& vars) -> ADdouble{
  	  return this->cost_func(vars) + this->alpha * this->penalty_func(vars);
   	  };
  	

  	 /*! @brief Solves the tree from different starting points
  	 
  	 @return an Rcpp::List with the results
    	*/
  	Rcpp::List solve_trees(void);
  

  	 /*! @brief Sets up mathematical moment
  	 PRoduces the OPtimFuns to send to the OptimHandler
  	 @note This function utilises the Visitor pattern. 	
  	    @param y the labels vector
  	    @param X_coeff the coefficients matrix  	
    	*/
  	OptimTraits::OptimFuns create_mathematical_model(const arma::vec& y, const arma::mat & X_coeff); 
  	
  	/*! @brief Sets up the optimiser member 
  	Sets variable sizes, constraints, variables and constraints bounds
  	@param y the labels vector
  	@param X_coeff the coefficients matrix  	
        */
  	void setup_optimiser(const arma::vec& y, const arma::mat & X_coeff);
  	/*! @brief Initialises, given a fixed seed, the variables vector
  	@param seed the random seed
  	@param vars the vector to alter.        
        */
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
