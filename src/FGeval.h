#include <cppad/ipopt/solve.hpp>
# include <cppad/utility/check_simple_vector.hpp>  // for compile-time checks
#include <vector>
#include <Rcpp.h>

using CppAD::AD;
// CheckSimpleVector<R, Vector>()
  
class FG_eval {
public:
  typedef CPPAD_TESTVECTOR( AD<double> ) ADvector; // works
  void operator()(ADvector& fg, const ADvector& x)
  {   
    assert( fg.size() == 3 );
    assert( x.size()  == 4 );
    
    // Fortran style indexing
    AD<double> x1 = x[0];
    AD<double> x3 = x[2];
    AD<double> x2 = x[1];
    AD<double> x4 = x[3];
    // f(x)
    fg[0] = x1 * x4 * (x1 + x2 + x3) + x3;
    // g_1 (x)
    fg[1] = x1 * x2 * x3 * x4;
    // g_2 (x)
    fg[2] = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4;
    //
    return;
  }
};