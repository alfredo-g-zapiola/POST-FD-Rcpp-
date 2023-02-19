#include "BasisObj.h"

FdHandler<BasisEnum::BSPLINE>::FdHandler(
    splines2::BSpline && bspline_basis):basis(std::move(bspline_basis)), N(21), quad(Simpson(), Mesh1D(Domain1D(0,1), N)),
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

double FdHandler<BasisEnum::BSPLINE>::operator()(const arma::mat& coef, 
                                               const unsigned i,
                                               const double t){
  this->basis.set_x(splines2::rvec{t});
  arma::mat basis_eval = this->basis.basis(true); // true to include intercept
  // (basis_eval * coefs_arma).eval()(0,0)
  return (basis_eval * coef.col(i)).eval()(0,0); 
  
};

arma::mat FdHandler<BasisEnum::BSPLINE>::compute_dissim_matrix(
    const arma::mat& X_coef){
  arma::mat dis_mat = arma::mat(X_coef.n_cols, X_coef.n_cols);
#ifndef MYNDEBUG
  std::cout << "Computing dissimilarity matrix" << std::endl;
#endif	
  
#if defined(PARALLELO) && defined(_OPENMP)
  int threads = omp_get_num_threads();
  Rcpp::Rcout << "open mp active. n threads recognised by default: " << threads << std::endl;
  Rcpp::Rcout << "note this package has 4 threads hardcoded if openmp enabled" << std::endl;
#pragma omp parallel for collapse(2) num_threads(1)
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
}; 


arma::mat FdHandler<BasisEnum::BSPLINE>::compute_features(
    const arma::mat & X_coef, unsigned n_feats){
#ifndef MYNDEBUG
  std::cout << "Computing features with n_feats = " << n_feats << std::endl;
#endif
  this->compute_basis_integrals(n_feats);
  
  // X_coef is df x n_samples ;
  // basis integrals is df x features
  // returns a matrix of n_samples * n_feats
  return  X_coef.t() * this->basis_integrals;
}


void FdHandler<BasisEnum::BSPLINE>::compute_basis_integrals(unsigned n_feats){
  
  this->basis_integrals.set_size(keep_bases_,n_feats);
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
    
#if defined(PARALLELO) and defined(_OPENMP)
#pragma omp parallel for shared(domain, quad, N) num_threads(1)
#endif
    for (unsigned k=0; k < keep_bases_; k++){
      auto meshgen = VariableSize(this->domain, weight_basis, this->N);
      // this->quad = Quadrature(Trapezoidal(), meshgen);
      basis_integrals(k, j) = quad.apply(
        [&weight_basis,k,this](double t){
          return weight_basis(t) * this->basis_function(k)(t);  // dot product
        }
      );
#ifndef MYNDEBUG
      std::cout << "Basis integral for basis " << k << ", feature " << j << \
        "is " << basis_integrals(k, j) << std::endl;
#endif
    }
  }
}
