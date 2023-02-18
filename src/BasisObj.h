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

enum class BasisEnum {  // TODO use for template, ALSO why class?
  BSPLINE = 0
};

template<BasisEnum b>
class FdHandler {

  virtual double operator()(const arma::mat& coef, unsigned i,double t) = 0;
  
  
};

template<>
class FdHandler<BasisEnum::BSPLINE>{
public:
  explicit FdHandler(splines2::BSpline && bspline_basis) : 
    basis(std::move(bspline_basis)), N(21), quad(Simpson(), Mesh1D(Domain1D(0,1), N)),
    keep_bases_(basis.get_spline_df()){
      #ifndef MYNDEBUG
      std::cout << "Constructing the basis handler" << std::endl;
      std::cout << "N nodes " << N << std::endl;
      std::cout << "number of bases " << keep_bases_ << std::endl;
      #endif
    
      arma::vec a_b = this->basis.get_boundary_knots();
      this->domain = Domain1D(a_b[0], a_b[1]);
      this->mesh = Mesh1D(domain, this->N);
      // this->basis_integrals.resize(keep_bases, keep_bases)
    
    
  };
  
  // TODO static polymorphism
  // evaluate at point t a sample 
  double operator()(const arma::mat& coef, unsigned i,double t){
    
    this->basis.set_x(splines2::rvec{t});
    arma::mat basis_eval = this->basis.basis(true); // true to include intercept
    // (basis_eval * coefs_arma).eval()(0,0)
    return (basis_eval * coef.col(i)).eval()(0,0); 
    
  };
  
  inline FunPoint basis_function(unsigned basis_idx){
    return [this, basis_idx](double t) -> double {
      this->basis.set_x(splines2::rvec{t});
      arma::mat basis_eval= this->basis.basis(true);
      return basis_eval.eval()(0, basis_idx);
      };
    };
  // TODO use weight function.  
  
  inline arma::mat compute_dissim_matrix(const arma::mat& X_coef){
    arma::mat dis_mat = arma::mat(X_coef.n_cols, X_coef.n_cols);
#ifndef MYNDEBUG
    std::cout << "Computing dissimilarity matrix" << std::endl;
#endif
    for (unsigned i = 0; i < X_coef.n_cols; i++){
#ifndef MYNDEBUG
      std::cout << "Column " << i << " of " << X_coef.n_cols<< std::endl;
#endif
      for (unsigned j = i+1; j < X_coef.n_cols; j++){
        dis_mat(i,j) = quad.apply([this, i,j, &X_coef](double t)-> double{
          return ((*this)(X_coef, i, t) - (*this)(X_coef, j, t) )*((*this)(X_coef, i, t) - (*this)(X_coef, j, t) );
          }
        );
        dis_mat(j,i) = dis_mat(i,j); 
#ifndef MYNDEBUG
        std::cout << "(i,j); " << i << " , "<< j << " is " << dis_mat(i,j)<< std::endl;
#endif
      }
      
    }
    return dis_mat;
  } 
  
  inline arma::mat compute_features(const arma::mat & X_coef, unsigned n_feats){
#ifndef MYNDEBUG
    std::cout << "Computing features with n_feats = " << n_feats << std::endl;
#endif
    this->compute_basis_integrals(n_feats);
    
    // X_coef is df x n_samples ;
    // basis integrals is df x features
    // returns a matrix of n_samples * n_feats
    return  X_coef.t() * this->basis_integrals;
  }
  
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
  void compute_basis_integrals(unsigned n_feats){
    this->basis_integrals.set_size(keep_bases_,n_feats);
    // TODO Pragma
    double I = domain.right() - domain.left();
    for (unsigned j = 0; j < n_feats; j++){
      double lb = this->domain.left() + j*I / (n_feats);
      double up = this->domain.left() + (j+1)*I/ (n_feats); 
#ifndef MYNDEBUG
      std::cout << "Current lb, up of integration for bases are: " << lb << " " << up << std::endl;
#endif
      auto weight_basis = [lb, up](double const&x) -> double{
        return (x < up and x>= lb) ? 1 : 0;  // step function
      };
      
      for (unsigned k=0; k < keep_bases_; k++){
        auto meshgen = VariableSize(this->domain, weight_basis, this->N);
        // this->quad = Quadrature(Trapezoidal(), meshgen);
        basis_integrals(k, j) = quad.apply(
          [&weight_basis,k,this](double t){
            return weight_basis(t) * this->basis_function(k)(t);
          }
        );
        #ifndef MYNDEBUG
        std::cout << "Basis integral for basis " << k << ", feature " << j << \
          "is " << basis_integrals(k, j) << std::endl;
          #endif
      }
    }
  }
    
  void compute_basis_squared_integrals(unsigned n_feats){
    
  }
  
  
};




#endif  // basis object