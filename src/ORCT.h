#ifndef ORCT_hh
#define ORCT_hh
// Note structure and var map depend on n_feats and depth

#include <functional>
#include <iostream> 
#include <vector>
#include <assert.h>   

#include "RcppArmadillo.h"
#include "helpers.h"
using Rcpp::_;  // named placeholder, used to create an Rcpp List


namespace fdpot{

/*! @brief Optimal Randomised Classification Tree
 * 
 * A class whose purpose is to deal with the ORCT's structure,
 * independently of whether it uses features of functionlal or multivariate data.
 * 
 */
struct ORCT{
  /*! @brief The type for the tree structure
   * 
   * Vector of a pair of vectors of unsigned:
   * It holds 
   * 
   */
  using StructureT=std::vector< std::pair<std::vector<unsigned>, std::vector<unsigned>> >;
  using TreeInfoT=std::vector<double>;
  
  
  /*! @brief Constructor
   * 
   * Given the input parameter, calculates other specs of the tree, namely
   * number of nodes, number of leaf noes and number of interior nodes
   * It also invokes the create_structure() method
   * @param depth_  the tree depth.
   * @param n_feats_ how many features used to perform a split at each node
   * @param n_labels_ how many different classes there are. 
   * @param gamma_  randomisation parameter
   */
  ORCT(const unsigned depth_, const unsigned n_feats_, const unsigned n_labels_,
       const double gamma_=512):
    depth(depth_), n_feats(n_feats_), n_labels(n_labels_), gamma(gamma_){
    this->n_nodes = helpers::n_nodes(depth);
    this->n_leaf_nodes = helpers::n_leaf_nodes(depth);

    this->n_int_nodes = n_nodes - n_leaf_nodes;
    #ifndef MYNDEBUG
      assert(this->n_int_nodes == helpers::n_nodes(depth-1) );
    #endif
    
    this->create_structure();
  }
  
  /*! @brief setup nodes and tree specs
   * 
   * This method fills the structure vector,
   * as well as the var_map
   * 
   */
  void create_structure(void); // 
  
  unsigned int depth;
  double gamma{512.};
  // note they could (or should) be made const but I woud have to rewrite the constructor
  // in the cpp file TODO
  unsigned n_nodes = 0, n_labels = 0, n_leaf_nodes = 0, n_int_nodes = 0, 
    n_feats = 0, n_samples = 0,n_vars = 0u;
  
  /*! @brief The tree structure
   * For its type, see the documentation this class: it is a vector of 
   * pairs of vectors of unsigned.
   * 
   * Each element of the outermost vector correspond to a node.
   * Within each node, there are two vectors of unsigned: the first one contains 
   * the left parent nodes' indices; the second one the right parent nodes' ones
   * 
   * @note Storing the parent nodes is crucial since the objective function for an 
   * ORCT works calculating the probability of falling in a certain node which
   * of course depends on the probability of falling in the nodes before such node
   */
  StructureT structure;
  
  /* @brief the bijection of variable 
   Maps each node index tau (unsigned) to the index of the first variable
   of such node in the Var
   @note it is a function pointer since it is completed in the create_structure()
   method
   @param
   @return
   */
  std::function<unsigned(unsigned)> var_map = nullptr;
  
  /* @brief computes the cumulative distribution function
   * 
   * Computes the cumulative distribution function value
   * @tparam VarT the variable type. In this library, two different types are used:
   * the ADdouble one (for CppAD compatibility for optimisation, see FdPot.cpp) and double for the 
   * predict method (see TODO)
   * 
   * @param x: the value at which we want to compute the cdf; in the ORCT it is 
   * the dot product of the features in a node and its variables plus the intercept
   *  @note See the report 
   */
  template<typename VarT>
  VarT cdf(const VarT x) const;

  
  /*! @brief computes the probability a statistical unit goes left
   * 
   * Computes the probability a statistical unit will go left, i.e. will go 
   * from current node tau to a node on the next level that is left of the current node
   * 
   * The features of the statistical unit (so in the functional case calculated ex ante)
   * are needed, as well as the whole variables vector. Hence this method utilises the
   * var_map (see var_map)
   * 
   * @tparam VarT the variable type. In this library, two different types are used:
   * the ADdouble one (for CppAD compatibility for optimisation, see FdPot.cpp) and double for the 
   * predict method 
   * enable_if could be used so that only these two types are used for template 
   * specialisation. For scalability it was not used.
   * 
   * @param feats the vector of the features for the statistical unit
   * @param vars the vector of all variables. Of
   * @param tau the node index (which node it is)
   * 
   */ 
  template<typename VarVecT> 
  typename VarVecT::value_type proba_go_left(const arma::rowvec & feats, const VarVecT &vars,
                         unsigned tau) const; 
    /*! @brief computes the probability a statistical falls on a given leaf
   * 
   * Computes the probability a statistical unit will fall on a given leaf,
   I.e. will follow a path of nodes that leads to such leaf
   * 
   * The features of the statistical unit (so in the functional case calculated ex ante)
   * are needed, as well as the whole variables vector. Hence this method utilises the
   * var_map (see var_map)
   * 
   * @tparam VarT the variable type. In this library, two different types are used:
   * the ADdouble one (for CppAD compatibility for optimisation, see FdPot.cpp) and double for the 
   * predict method 
   * enable_if could be used so that only these two types are used for template 
   * specialisation. For scalability it was not used.
   * 
   * @param feats the vector of the features for the statistical unit
   * @param vars the vector of all variables. Of
   * @param tau the leaf node index
   * 
   */ 
  template<typename VarVecT>
  typename VarVecT::value_type proba_fall_leaf(const arma::rowvec&feats, const VarVecT & vars,
                           unsigned tau) const;
  /*! @brief predict the labels (both probability and actual value)
  
  @param feats the features computed from the sample
  @param the fitted variables after optimisation
  @return an Rcpp list with probabilities of belonging to the labels and the labels with maximum probability
  */
  Rcpp::List predict(const arma::mat& coefs, const arma::vec& vars) const;

  
  
};

// clarify why I did not use constexpr template for exp.

template<typename VarVecT>
typename VarVecT::value_type ORCT::proba_go_left(const arma::rowvec & feats,  
                             const VarVecT & vars,unsigned tau) const{
  using VarT = typename VarVecT::value_type;
  VarT val = 0.;
      
  unsigned var_idx = this->var_map(tau); 
  
  #ifndef MYNDEBUG
    // assert tau does not belong to last level
    assert(tau < this->n_int_nodes);
    // assert the var index is within the indices for the 
    assert(var_idx < this->n_int_nodes * (this->n_feats + 1));
    // assert N-elements matche sn_feats
    assert(feats.n_elem == this->n_feats);
    // std::cout << "feats.n_elem " << feats.n_elem << std::endl;
    // std::cout << "N_feats  " << this->n_feats << std::endl;
    
  #endif
      
  
  for (unsigned i = 0; i < feats.n_elem; i++){  // perform dot product
    val += feats(i) * (vars[var_idx++]);  // cannot use iterators with pragma (dependent)
  }
  // normalise by number of features
  val /= feats.n_elem; 
  #ifndef MYNDEBUG
    assert(var_idx == this->var_map(tau) + this->n_feats);
  #endif
  // subtract the mu variable (intercept)
  val -= vars[var_idx];
  
  return cdf<VarT>(val);
  
}


/*! @brief Probability the SU falls on a leaf

@tparam VarVecT the type of the variable vector, in this project either arma::vec or OptimTraits::ADvector


*/
template<typename VarVecT>
typename VarVecT::value_type ORCT::proba_fall_leaf(const arma::rowvec& feats,
                                          const VarVecT & vars,
                                          unsigned tau) const{
  using VarT = typename VarVecT::value_type;
  
#ifndef NDEBUG
  // assert that tau, the node index corresponds to one of the leaves and not
  // an interior node
  assert(tau >= helpers::n_nodes(depth - 1) and tau < this->n_nodes );
#endif 
  
  auto& leaf = structure[tau];
  VarT proba = 1.;
  
  for (auto l: leaf.first)  // leaf.first is array of parents of nodes that went left
    proba *= proba_go_left<VarVecT>(feats, vars, l);
  
  for (auto r: leaf.second) // leaf.first is array of parents that went right
    proba *= (1 - proba_go_left<VarVecT>(feats, vars, r));
  
  return proba;  
  
  
}

}  // namespace fdpot

#endif // of the ORCT header file
