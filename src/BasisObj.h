#ifndef BASIS_OBJ_HH
#define BASIS_OBJ_HH

#include <splines2Armadillo.h>

// for integration
#include "numerical_integration.hpp"
#include "Adams_rule.hpp"


using namespace apsc::NumericalIntegration;
using apsc::NumericalIntegration::FunPoint;
using apsc::NumericalIntegration::Simpson;
using namespace Geometry;

/*! @brief Enumeration for the different bases
 * 
 * An enum class to keep track of the possible bases the FdHandler template
 * may specialise in
 * @note as explained in the report, only the BSPLINE basis type is currently developed. 
 */
enum class BasisEnum {
  BSPLINE = 0
};

/*! @brief Handle Functional Data feature and dissimilarity
 * 
 * The base template method for the Functional Datum.
 * 
 * @tparam b the Basis type (given by BasisEnum) the class can specialise int.
 * @note as explained in the report, only the BSPLINE basis type is currently developed. 
 */
template<BasisEnum b>
class FdHandler {
/*! @brief evaluate ith functioanl datum at point t
 * 
 
 @note It is not const on purpose, since a workaround to use the splines2::BSpline
 object was to use the non-const set_x method.
 @param coef: the matrix of coefficients fitted in the smoothing
 @param i: the index of the statistical unit (which function)
 @param t: point at which the function should be evaluated
 @return the evaluation of the function, of course a double.
 */
  virtual double operator()(const arma::mat& coef, const unsigned i,const double t) = 0;
  
  
};

template<>
class FdHandler<BasisEnum::BSPLINE>{
public:
/*! @brief Constructor

@param bspline_basis an rvalue, since the handler receives it from outside but then owns it.

*/
  explicit FdHandler(splines2::BSpline && bspline_basis);
  
	 
 /*! @brief evaluate ith functioanl datum at point t
 * 
 This function calls the basis_function method to obtain the evaluation of the bases at pint t.
 then performs the dot product with the coefficients to obtain the evaluation.
 
 @note It is not const on purpose, since a workaround to use the splines2::BSpline
 object was to use the non-const set_x method.
 @param coef: the matrix of coefficients fitted in the smoothing
 @param i: the index of the statistical unit (which function)
 @param t: point at which the function should be evaluated
 @return the evaluation of the function, of course a double.
 */
  double operator()(const arma::mat& coef, const unsigned i, const double t);
/*! @brief Obtain the callable of a basis function
   evaluates the basis contained in the splines2::BSpline member at the given point
   @param basis_idx the index of the basis
   @return a function that when called returns the value of such basis at the passed point
 */
  inline FunPoint basis_function(unsigned basis_idx){
    return [this, basis_idx](double t) -> double {
      this->basis.set_x(splines2::rvec{t});
      arma::mat basis_eval= this->basis.basis(true);
      return basis_eval.eval()(0, basis_idx);
      };
    };
  /*! @brief Compute dissimilarity matrix of funcitonal data
   Calls the operator() to 
   Then applies the dot product
   @param X_coef the 
   @return a function that when called returns the value of such basis at the passed point
 */
  arma::mat compute_dissim_matrix(const arma::mat& X_coef);
  
 /*! @brief Compute the features from functional data
   Performs the dot product of the functions with step functions.
   Once the integrals of the bases * step functions are obtained, the dot product of that value
   and the coefficients is performed (see report).
   
   @param X_coef the coefficient matrix
   @param n_feats how many features (if 2, 1 step function from 0 to 0.5 and another from 0.5 to 1)
   @return a function that when called returns the value of such basis at the passed point
 */
  arma::mat compute_features(const arma::mat & X_coef, unsigned n_feats);
  
private:
  // members
  double a, b;
  unsigned n_feats=0;
  splines2::BSpline basis;
  const unsigned N;
  Domain1D domain;
  Mesh1D   mesh;
  Quadrature quad;
  unsigned keep_bases_;
  arma::mat basis_integrals;
  void compute_basis_integrals(unsigned n_feats);
    
  // void compute_basis_squared_integrals(unsigned n_feats){}; // TODO 
  
  
};




#endif  // basis object
