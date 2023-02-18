#ifndef OPTIM_HANDLER_HH
#define OPTIM_HANDLER_HH
#include <cppad/ipopt/solve.hpp>
#include <functional>
#include <memory>
#include <vector>

namespace fdpot{
class FdPot;  // forward declaration

using CppAD::AD;
using ADdouble = AD<double>;

/*!
  @brief Keeps the types for optimisation
 
 @note See https://www.coin-or.org/CppAD/Doc/ipopt_solve.htm
 */
struct OptimTraits{

/**  ADvector: the type of vector supported by CppAD */
 using ADvector=std::vector<AD<double>>; //using Dvector= ADvector;
 
  /**  OptimFuns: the type to hold optimisation functions. wrapper for function pointer*/
 using OptimFuns=std::function<void(ADvector& fg, const ADvector& x)>;
 
 /**  Dvector: the vector that will hold the variables changed during optimisation */
 using Dvector=arma::vec;  // could also be std::vector

  
};

/*! @brief An abstract class to evaluate function and gradients

 @description CppAD has a very simple (and ugly) example for an interface to Ipopt.
 I decided to make an abstract class so that the functor can be set once the tree is built
 so that the objective function and the variables vector (which depends on the depth, number of feats and number of labels)
can be passed later.
  @note an abstract class of the one seen here: https://www.coin-or.org/CppAD/Doc/ipopt_solve_get_started.cpp.htm
*/
class Abstract_FG_eval: public OptimTraits {
public:
  //  typedef CPPAD_TESTVECTOR( AD<double> ) ADvector;
  
  /*! @brief typedef as required (see the link)
  
  @note see  https://www.coin-or.org/CppAD/Doc/simplevector.htm
   */
  typedef std::vector<AD<double>> ADvector; 
  
  /*! @brief 
  @param fg vector with objective function and constraint function values
  @param x  vector with the optimising varaibles on which fg depends
  */
  virtual void operator()(ADvector& fg, const ADvector& x) = 0; // make it abstract
};

struct FG_eval: public Abstract_FG_eval{
  /*! @brief 
   * @param optimfuns_ the optimisation functions 
   */
   FG_eval(const OptimFuns& optimfuns_): opt_funs{optimfuns_}{};
   FG_eval(void) = default;
  
  /*! @brief 
   @description Evaluates the objective function and constraint functions utilising a private member that contains the optimization functions
   @param fg vector with objective function and constraint function values
   @param x  vector with the optimising varaibles on which fg depends
   */
  void operator()(ADvector& fg, const ADvector& x) override final {
    this->opt_funs(fg, x); // 
    return; 
  };
  
  OptimFuns opt_funs = nullptr;
  
};

/*! @brief Interface class for optimisation
 * 
 * @description This class provides the methods to perform an optimisation, 
 * namely for variable creation, objective function and constraints definition 
 * and solving
 
*/
struct OptimHandler: public OptimTraits{
public:
  //using Dvector = std::vector<double>;
  
  /*! @brief Vector type for variables and constraints functions values
   */
  using Dvector = arma::vec;
  /*! @brief Default void constructor
   
   Takes no arguments since the optimisation procedure occurs calling the
   methods of this class. Indeed, only once the .fit() method of the FdPot is
   called can the optimiser be set up (number of variables, ecc.)
   * 
   */
  OptimHandler(void) = default;  
  
  /*! @brief Number of variables in the opt problem
   * @note Updated when calling the public methods of this class
   */
  unsigned n_vars = 0u;
  
  /*! @brief Number of constraints in the opt problem
   * @note Updated when calling the public methods of this class
   */
  unsigned n_constraints = 0u;
  
  /*! @brief "Ugly interface for Ipopt"
   * @note The required class by Ipopt and CppAD to solve an opt. problem this
   * class encapsulates to provide a better interface
   */
  FG_eval  fg_eval= FG_eval();  // std::unique_ptr<FG_eval>
  
  /*! @brief Vector holding the optimising variables
   The vector holding the variables that are updated during the optimisation
   */
  Dvector variables;
  
  
  /*! @brief bounds for vars and constraints
    xl: lower bounds for the variables (first element matches the first variable)
    xu: upper bounds " ""
    gl: lower bounds for the constraints functions
    gu upper bounds for the constraints function
  */
  Dvector  xl, xu, gl, gu; 
  
  /*! @brief member that stores the solution
  Once the optimisation is done, the solution will be stored in this member.
  See the solve method
  */
  CppAD::ipopt::solve_result<Dvector> solution;
  
  /*! @brief internal adjustments given number of vars and constrs
   * 
   * Given the number of variables and of constraints, updates the members
   * that store these values; also resizes the variables vector
   */
  inline void create_variables(const unsigned n_vars_,
                               const unsigned n_constrs_){
    
    this->n_vars = n_vars_;
    this->n_constraints = n_constrs_;
    variables.resize(n_vars); 

    
  }
  /*! @brief Update all members that deal with optimisation
   * All the math program specs should be specified here: 
   * optimisation functions (objective and constraints) as well as the 
   * bounds for variables and constraint functions
   */
  inline void set_math_program(const OptimFuns optimfuns,
                               Dvector&& xl,
                               Dvector&& xu,
                               Dvector&& gl,
                               Dvector&& gu){
    this->fg_eval = FG_eval(optimfuns);
    this->xl  = std::move(xl);
    this->xu = std::move(xu);
    this->gl= std::move(gl);
    this->gu= std::move(gu);
    // TODO options
  }
  /*! @brief perform the optimisation
   
   This method utilises the FG_eval class, aggregated in an unique pointer 
   in this class
   * @note options are left as the default ones in the CppAD example,
   * except for the numeric tolerance
   */
  inline void solve(void){ // TODO assert sizes
    // solve the problem
    std::string options;
    // options += "Integer print_level  0\n";  // disable verbose
    options += "String  sb           yes\n";
    // maximum number of iterations
    // options += "Integer max_iter     10\n";
    // approximate accuracy in first order necessary conditions;
    // see Mathematical Programming, Volume 106, Number 1,
    // Pages 25-57, Equation (6)
    options += "Numeric tol          1e-12\n";
    // derivative testing
    options += "String  derivative_test            second-order\n";
    // maximum amount of random pertubation; e.g.,
    // when evaluation finite diff
    options += "Numeric point_perturbation_radius  0.\n";
#ifndef MYNDEBUG
    std::cout << "Variables: " << variables << std::endl;
#endif
#ifdef DEV
    Rcpp::Rcout << "Variables: " << variables << std::endl;
#endif
        CppAD::ipopt::solve<Dvector, FG_eval>(options, variables, xl, 
                                                      xu, gl, gu, 
                                                      this->fg_eval, solution);
  }; 
  
};
  
  
    
  
} // namespace fdpot


#endif 