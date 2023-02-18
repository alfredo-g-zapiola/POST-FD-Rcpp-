#include "ORCT.h"
namespace fdpot{

template<>
double ORCT::cdf<double>(const double x) const {
  return 1 / (1 + std::exp(-x * this->gamma));
}

void ORCT::create_structure(void){
  // map for the left parent and right parent nodes
  // bijection for variables 
  this->structure.resize(this->n_nodes);
  
#ifndef MYNDEBUG
  std::cout << "Structure size is: " << structure.size() << std::endl;
#endif 
  
  //  structure.push_back(std::make_pair(std::vector<unsi))
  structure[0] = std::make_pair(std::vector<unsigned>(), std::vector<unsigned>());  
  // the current node we will edit
  unsigned node = 1;
  
  for (unsigned level = 1; level <= depth; level++){
    // first node index of previous level coincides with the number of nodes
    // of its anterior level's number of nodes
    unsigned prev_idx  = helpers::n_nodes(level - 2);
    // expand from each node of previous level
    for (unsigned j = prev_idx; j < helpers::n_nodes(level - 1); j++){
      // extract ancestors of left and right nodes
      auto left_nodes = std::vector<unsigned>(structure[j].first.cbegin(), 
                                              structure[j].first.cend());
      auto right_nodes = std::vector<unsigned>(structure[j].second.cbegin(),
                                               structure[j].second.cend());
      // TODO initialise with zeroes to allocate space
      structure[node] = std::make_pair(left_nodes, right_nodes);
      // since we went to the left, j is a left parent node, so update first
      structure[node].first.push_back(j);  
      std::cout << "Node " << node  << std::endl;
      for (const auto v: structure[node].first)
        std::cout << "left parent value "  << v << std::endl;
      for (const auto v: structure[node].second)
        std::cout << "right parent value "  << v << std::endl;
      node++;

      // go right
      // now we use move since we won't need the vecs anymore, and we save 
      // some computations
      structure[node] = std::make_pair(std::move(left_nodes), std::move(right_nodes));
      // we went to the right, update second
      structure[node].second.push_back(j);
      std::cout << "Node " << node  << std::endl;
      for (const auto v: structure[node].first)
        std::cout << "left parent value "  << v << std::endl;
      for (const auto v: structure[node].second)
        std::cout << "right parent value "  << v << std::endl;
      node++;
      
      // Now we update the var_map: for each node tau maps the index of its first
      // variable
      this->var_map = [this] (const unsigned tau) -> unsigned{
        if (tau == 0)
          return 0;
        else if (tau < this->n_int_nodes)
          return (this->n_feats + 1) * (tau);  // number of features + intercept
        else{
          // compute index of last 
          unsigned cum_int_idx = (this->n_int_nodes * (this->n_feats + 1));
          return  cum_int_idx + this->n_labels * (tau - this->n_int_nodes);
        }
      };
    }
  }
  
  this->n_vars = (n_feats + 1) * n_int_nodes + n_labels * n_leaf_nodes;
};


Rcpp::List ORCT::predict(const arma::mat& feats, 
                         const arma::vec& vars) const{
  arma::mat probs_mat(feats.n_rows, this->n_labels);
  arma::vec labels_vec(feats.n_rows);
  
  for (unsigned i = 0; i < feats.n_rows; i++){  // for each new statistical unit
    
    for(unsigned k = 0; k < this->n_labels; k++){  // for each label
      
      double prob_k = 0.;  // initialise the probability of kth label
      // now iterate are leaves
      for (unsigned leaf = this->n_int_nodes;  // the 
                    leaf < this->n_nodes; leaf++){
        // first_class var in this leaf index
        unsigned c_kt = this->var_map(leaf) + k;
        // add the probability of falling on current leaf times \
        // probability kth class is chosen
        prob_k += this->proba_fall_leaf(feats.row(i), vars, leaf) * vars[c_kt];
      }
      
      probs_mat(i, k) = prob_k;
    }  // end for k
   
    #ifndef MYNDEBUGs
        double EPS{.001};
        // validate if the sum of probabilities is 1
        assert( std::abs(arma::sum( probs_mat.row(i) ) - 1) < EPS);
    #endif
    
    unsigned best_label = probs_mat.row(i).index_max();
    labels_vec(i) = best_label;
    
  } // end for i
  return Rcpp::List::create(_("predicted_labels_probs") = probs_mat,
                            _("predicted_labels") = labels_vec);
  };
}; // namespace fdpo

