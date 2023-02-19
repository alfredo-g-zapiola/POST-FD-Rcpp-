#ifndef FDPOT_HELPERS
#define FDPOT_HELPERS
#include <cmath>

namespace helpers{	
/**
 * Given a level, obtains the cumulative sum of nodes present in the tree
 * 
 * @param level the depth of the tree
 * @return number of total nodes are already present when increasing dept
 */
inline unsigned cum_tree_sum(unsigned level){
  if (level == 0)
    return 0;
  else{
    unsigned sum = 0;
    for (unsigned i = 0; i <= level; i++)
      sum += std::pow(2, i);
    return sum;
  }
};

/**
 * Given a level, obtains the cumulative sum of nodes until that level (excluded)
 * 
 * @param level of the tree, smaller or equal than the depth
 * @return number of total nodes are already present when increasing dept
 */
inline unsigned n_nodes_until_level(unsigned level){
  unsigned sum = 0;
  for (unsigned i = 0; i < level; i++)
    sum += std::pow(2, i);
  return sum;
}

/**
 * Obtain the total number of nodes given the tree depth
 *
 * @param depth of the tee
 * @return total number of nodes
 */
inline unsigned n_nodes(unsigned depth){
  return std::pow(2, depth + 1) -1 ;
}

/**
 * Obtain the number of leaf nodes nodes given the tree depth
 *
 * @param depth of the tee
 * @return total number of nodes
 */
inline unsigned n_leaf_nodes(unsigned depth){
  if (depth == 0)
    return 1;
  else{
    return n_nodes(depth) - n_nodes(depth-1);
  }
  
}

} //namespace helpers

#endif 
