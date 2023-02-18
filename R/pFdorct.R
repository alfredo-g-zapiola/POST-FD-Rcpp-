#' Build and fit an FD-POT that classifies functional data
#'@description Given a vector of labels and an fdSmooth object (smoothed functional data), builds and fits a penalised optimal randomised decision tree
#'
#'
#'@param y vector of integers with the labels. They MUST be numbered starting from 0
#'@param X the smoothed functional dats object of class fdSmooth
#'@param depth the classificatin tree depth
#'@param alpha the weight given in the objective function (stronger alpha, higher penalty for dissimilarity in each leaf node)
#'@param similarity.method the method used to obtain the dissimilarity between two functional data. Currently only L2 norm of the difference supported
#'@param n_feats the number of features to compute (integrals) for each functional datum
#'@param n.solve how many optimisations to carry out (each from a different starting point)
#'@param m_cost the misclassification cost (best left at default value)
#'@param gamma the randomisation factor (best left at default value)
#'@param seed the random seed
pFdorct <- function(y, X, basis.degree, depth = 2, alpha = .5, similarity.method="d0.L2", 
                    n_feats=10, n.solve = 20,gamma=512, seed=21071865){
  # TODO ask parameters for degree
  if (! class(X) == "fdSmooth"){
    stop("X must be of fdSmooth class")
  }
  basis.type <- X$fd$basis$type
  if (basis.type == "bspline"){
    print("Calling the rcpp function")
    # call the Rcpp function which instantiates the C++ class and fits the model
    tree.datalist <- pFdorct_Rcpp(y, X$fd$coefs, X$argvals, as.integer(X$df),
                              basis.degree,
                              "BSpline",
                              depth=depth,
                              alpha=alpha,
                              similarity_method=similarity.method,
                              n_feats = n_feats,
                              n_solve=n.solve,
                              gamma=gamma,
                              seed=seed) 
  }
  else{
    stop("only the bspline basis type is currently supported")
  }
  # prepare the info dispatch
  res = tree.datalist
  # make it an S3 class
  class(res) = "p.fdorct"
  return(res)
}

#'
#'@description
#'
#'@param
#'
predict.pFdorct <- function(model, X_fd_new){
  if (! class(model) == "p.fdorct")
    stop("model must be of p.fdorct class")
  if (! class(X_fd_new) == "fdSmooth"){
    stop("X must be of fdSmooth class")
  }
  if (X_fd_new$fd$basis$type != "bspline")
    stop("only the bspline basis type is supported")
  
  return(predict_FdPot_Rcpp(model, X_fd_new$fd$coefs))
  
}

#summary.pFdorct <- function(model, type = "class"){}
  
