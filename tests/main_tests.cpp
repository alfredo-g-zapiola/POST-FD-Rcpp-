#include <iostream>
# include <cppad/ipopt/solve.hpp>
#include "FdPot.h"
#include "chrono.hpp"

//////////////////////
// TEST 1
//////////////////////

// This is the example found on: https://www.coin-or.org/CppAD/Doc/ipopt_solve_get_started.cpp.htm
// I perform this test because I want to test if I can avoid doing the references
// (see below) to the elements of ADvector x when calling FG_eval


namespace {
using CppAD::AD;
  
  class FG_eval {
  public:
    //  typedef CPPAD_TESTVECTOR( AD<double> ) ADvector;
    // let's try with another container:
    // note: https://www.coin-or.org/CppAD/Doc/simplevector.htm
    typedef std::vector<AD<double>> ADvector; // Works too!
    
    // using ADvector=Rcpp::NumericVector; // as said by Cppad, any SimpleVector type works
    
    void operator()(ADvector& fg, const ADvector& x)
    {   assert( fg.size() == 3 );
      assert( x.size()  == 4 );
      // in the example, references to the vars were made here.
      // instead, I access directly the elements of the vector
      // f(x)
      fg[0] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
      // g_1 (x)
      fg[1] = x[0] * x[1]* x[2] * x[3];
      // g_2 (x)
      fg[2] = x[0] * x[0] + x[1] *x[1] + x[2] *x[2] + x[3] * x[3];
      //
      return;
    }
  };
}

bool get_started(void)
{   bool ok = true;
  size_t i;
  // typedef CPPAD_TESTVECTOR( double ) Dvector;
  typedef std::vector<double> Dvector;
  
  // number of independent variables (domain dimension for f and g)
  size_t nx = 4;
  // number of constraints (range dimension for g)
  size_t ng = 2;
  // initial value of the independent variables
  Dvector xi(nx);
  xi[0] = 1.0;
  xi[1] = 5.0;
  xi[2] = 5.0;
  xi[3] = 1.0;
  // lower and upper limits for x
  Dvector xl(nx), xu(nx);
  for(i = 0; i < nx; i++)
  {   xl[i] = 1.0;
    xu[i] = 5.0;
  }
  // lower and upper limits for g
  Dvector gl(ng), gu(ng);
  gl[0] = 25.0;     gu[0] = 1.0e19;
  gl[1] = 40.0;     gu[1] = 40.0;
  
  // object that computes objective and constraints
  FG_eval fg_eval;
  
  // options
  std::string options;
  // turn off any printing
  options += "Integer print_level  0\n";
  options += "String  sb           yes\n";
  // maximum number of iterations
  options += "Integer max_iter     10\n";
  // approximate accuracy in first order necessary conditions;
  // see Mathematical Programming, Volume 106, Number 1,
  // Pages 25-57, Equation (6)
  options += "Numeric tol          1e-6\n";
  // derivative testing
  options += "String  derivative_test            second-order\n";
  // maximum amount of random pertubation; e.g.,
  // when evaluation finite diff
  options += "Numeric point_perturbation_radius  0.\n";
  
  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;
  
  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, xi, xl, xu, gl, gu, fg_eval, solution
  );
  //
  // Check some of the solution values
  //
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  //
  double check_x[]  = { 1.000000, 4.743000, 3.82115, 1.379408 };
  double check_zl[] = { 1.087871, 0.,       0.,      0.       };
  double check_zu[] = { 0.,       0.,       0.,      0.       };
  double rel_tol    = 1e-6;  // relative tolerance
  double abs_tol    = 1e-6;  // absolute tolerance
  for(i = 0; i < nx; i++)
  {   ok &= CppAD::NearEqual(
    check_x[i],  solution.x[i],   rel_tol, abs_tol
  );
    ok &= CppAD::NearEqual(
      check_zl[i], solution.zl[i], rel_tol, abs_tol
    );
    ok &= CppAD::NearEqual(
      check_zu[i], solution.zu[i], rel_tol, abs_tol
    );
  }
  
  return ok;
}


//////////////////////
// TEST 2
//////////////////////

void test_integrals(void){
  auto X_argvals = arma::linspace(0, 1, 365);
  arma::vec boundary_knots({0,1});
  unsigned df{20u}, degree{3u};
  
  splines2::BSpline bs_obj = splines2::BSpline(X_argvals, df, degree,
                                               boundary_knots);
  FdHandler<BasisEnum::BSPLINE> fd_handler(std::move(bs_obj));
  // Now generate random data
  std::cout << "Generating test matrix C" << std::endl;
  arma::mat basis_coefs = arma::mat(df, 10).ones();
  std::cout << "Basis coefs: i.e. matrix C" <<  basis_coefs << std::endl;
  // basis_coefs.print();
  std::cout << "Computing features:" << std::endl;
  auto ints = fd_handler.compute_features(basis_coefs, 3);
  std::cout << ints << std::endl;
  std::cout << "Computing dissim matrix" << std::endl;
  auto disim_mat = fd_handler.compute_dissim_matrix(basis_coefs);
  std::cout << disim_mat << std::endl;
  
  //auto featmat = fd_handler.compute_features(10);
  // std::cout << featmat << std::endl;
  
  std::cout << std::endl;
  std::cout << "------------------ End of test 2 ------------------" << std::endl;
  
}


//////////////////////
// TEST 3
//////////////////////
void test_3(void){
  std::cout << "Obtaining basis object" << std::endl;
  arma::vec X_argvals = arma::linspace(0, 1, 100);
  std::cout << "ARGAVLs: " << X_argvals[0] << " and " << X_argvals[X_argvals.size()-1] << std::endl;
  arma::vec boundary_knots{X_argvals[0], X_argvals[X_argvals.size()-1]};
  std::cout << "Printing boundary knots: \n" << boundary_knots << std::endl;
  const unsigned DF = 20u;
  auto basis = splines2::BSpline(X_argvals, DF, 3u, boundary_knots);
  std::cout << "Instantiating tree" << std::endl;
  
  fdpot::FdPot tree = fdpot::FdPot(std::move(basis), 2, 10, 12, 2u, .1, 41703192, 512);
  std::cout << "Finished building tree" << std::endl;
  arma::vec y = arma::vec(10);
  for (unsigned i = 0; i < y.size(); i++){
    y(i) = i % 2 == 0 ? 1 : 0;
  }
  
  arma::mat test_coef_matrix;
  std::cout << "Loading from csv" << std::endl;
  test_coef_matrix.load("./test_coef_matrix.csv", arma::csv_ascii);
  std::cout << "Printing" << std::endl;
  std::cout << "Printing 0,0 " << test_coef_matrix(0,1) << std::endl;
  for (int i = 0; i < test_coef_matrix.n_rows; i++)
    std::cout << i << ", 2): " << test_coef_matrix(i, 2) << std::endl;
  std::cout << "fITTING TREE WITH " << y << std::endl;
  auto result = tree.fit(y, test_coef_matrix);
  // auto result = tree.fit(y, arma::mat(20, 10).ones() );
  std::cout << "Finished fitting tree" << std::endl;
  std::cout << "------------- END OF TEST 3 --------------" << std::endl;
  
  
}

void time_integrals(void){
  Timings::Chrono myclock;
  
  // setup the splines
  auto X_argvals = arma::linspace(0, 1, 365);
  arma::vec boundary_knots({0,1});
  unsigned df{20u}, degree{3u};
  
  splines2::BSpline bs_obj = splines2::BSpline(X_argvals, df, degree,
                                               boundary_knots);
  //auto basis = splines2::BSpline(X_argvals, X_basis_df, X_basis_degree,
  //                       boundary_knots);
  arma::mat basis_coefs = arma::mat(df, 10).ones();
  FdHandler<BasisEnum::BSPLINE> fd_handler(std::move(bs_obj));
  
  auto ints = fd_handler.compute_features(basis_coefs, 3);
  
  myclock.start();
  for (unsigned i = 0; i < 100; i++)
    auto disim_mat = fd_handler.compute_dissim_matrix(basis_coefs);
  myclock.stop();
  std::cout << "Timing of the dissimilarity matrix: " << myclock << std::endl;
  
  // Now measure the time for
  arma::mat test_coef_matrix;
  test_coef_matrix.load("./test_coef_matrix.csv", arma::csv_ascii);
  Timings::Chrono myclock2;
  myclock.start();
    for (unsigned i = 0; i < 100; i++)
      fd_handler.compute_features(test_coef_matrix, 4);
  myclock.stop();
  std::cout << "Timing of the features computation: " << myclock << std::endl;
  
}

int main(void){
  std::cout << "Test 1: Cpp interface" <<  std::endl;
 bool state = get_started();
  std::cout << "State is (should be 1): " << state << std::endl;
  
  std::cout << "Test 2: integrals" << state << std::endl;
  test_integrals(); 
  std::cout << "Test 3: timing integrals" << std::endl;
  time_integrals();
std::cout << "Test 4: tree" << state << std::endl;
 test_3();

  
#ifdef _OPENMP
  int threads = omp_get_num_threads();
  std::cout << "open mp active. n threads " << threads << std::endl;
#endif
  
}