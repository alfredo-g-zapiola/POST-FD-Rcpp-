library(FdPot)

df.X <- read.csv("data/X_canada.csv", header = F)

y <- read.csv("data/y_canada.csv", header=F)
train.idx <- as.matrix(read.csv("data/train_indices.csv", header=F))
X.train <- df.X[train.idx,]
df.X <- t(df.X)
X.train <- t(X.train)
y.train <- y[train.idx]

m <- 5           # spline order 
degree <- m-1    # spline degree 
nbasis = 20
basis <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis, norder=m)
time = seq(0, 1, length.out = 365)
Xsp <- smooth.basis(argvals=time, y=X.train, fdParobj=basis)
# save to file  
write.table(matrix(Xsp$fd$coefs, ncol=dim(Xsp$fd$coefs)[2],
                 nrow=(dim(Xsp$fd$coefs)[1])),
          file = "data/test_coef_matrix.csv", sep=",",
          row.names=F,
          col.names=F)


# ocnstruct with wrong number of labels
tryCatch(
  pFdorct(sample(c(0,1), replace=TRUE, size=10), Xsp, degree),
  error = function(e){print(e)}
)
# this throws an error 
n.labels <- dim(Xsp$fd$coefs)[2]
set.seed(20)
a.well <-  pFdorct(y.train, Xsp, degree, depth = 2, alpha = .1, n.solve = 5,
                   n_feats = 4)
print(paste0("Printing objective function values", 
             a.well$fit_results$obj_func_vals))
print("Print cost function value")
print(a.well$fit_results$cost_func_vals)
predis <- predict_FdPot_Rcpp(a.well, Xsp$fd$coefs, a.well$fit_results$best_tree_idx)
predis
preds = predict.pFdorct(a.well, Xsp)
sum(preds$predicted_labels_probs[1,])
vars = a.well$fit_results$best_variables
mean(preds$predicted_labels == y.train)  # accuracy


