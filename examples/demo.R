library(FdPot)
# setwd("path_to_fdpot")
# see documentation
?`FdPot-package`
?FdPot::pFdorct_Rcpp
?FdPot::pFdorct

# 1. load data
df.X <- read.csv("data/X_canada.csv", header = F)

y <- read.csv("data/y_canada.csv", header=F)
train.idx <- as.matrix(read.csv("data/train_indices.csv", header=F))
X.train <- df.X[train.idx,]
df.X <- t(df.X)
X.train <- t(X.train)
y.train <- y[train.idx]

# SETUP spline bassis
m <- 5           # spline order 
degree <- m-1    # spline degree 
nbasis = 20
basis <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis, norder=m)
time = seq(0, 1, length.out = 365)

# smooth functional data
Xsp <- smooth.basis(argvals=time, y=X.train, fdParobj=basis)
x11()
plot.fd(Xsp$fd)
graphics.off()

# now fit the tree
tree <-  pFdorct(y.train, Xsp, degree, depth = 2, alpha = .01, n.solve = 5,
                   n_feats = 4) 
print("Printing obj function values")
print(tree$fit_results$obj_func_vals)
print("Print cost function values")
print(tree$fit_results$cost_func_vals)
print("printing variables")
vars = a.well$fit_results$best_variables

# perform prediction:
preds = predict.pFdorct(tree, Xsp, tree$fit_results$best_tree_idx)
sum(preds$predicted_labels_probs[1,])

print("accuracy train set:")
print(mean(preds$predicted_labels == y.train))  # accuracy train set

# now fit the tree WITHOUT penalty
tree <-  pFdorct(y.train, Xsp, degree, depth = 2, alpha = 0, n.solve = 5,
                 n_feats = 4) 
print("Printing obj function values")
print(tree$fit_results$obj_func_vals)
print("Print cost function values")
print(tree$fit_results$cost_func_vals)
print("Print penalty function values")
print(tree$fit_results$cost_func_vals)
print("printing variables")
vars = a.well$fit_results$best_variables

# perform prediction:
preds = predict.pFdorct(tree, Xsp, tree$fit_results$best_tree_idx)

print("accuracy train set:")
print(mean(preds$predicted_labels == y.train))  # accuracy train set

