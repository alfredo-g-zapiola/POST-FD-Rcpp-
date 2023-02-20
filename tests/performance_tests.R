library(FdPot)
library(splines2)
library(microbenchmark)

# test 1: efficiency of the integrals
#====
# read data
df.X <- read.csv("data/X_canada.csv", header = F)
y <- read.csv("data/y_canada.csv", header=F)
train.idx <- as.matrix(read.csv("data/train_indices.csv", header=F))
X.train <- df.X[train.idx,]
df.X <- t(df.X)
X.train <- t(X.train)
y.train <- y[train.idx]

# set the functional data
m <- 5           # spline order 
degree <- m-1    # spline degree 
nbasis = 20
basis <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis, norder=m)
time = seq(0, 1, length.out = 365)
Xsp <- smooth.basis(argvals=time, y=X.train, fdParobj=basis)

# extract information on the func data
coefs = Xsp$fd$coefs[,1]
internal_knots = FdPot::get_bspline_internal_knots(coefs, Xsp$argvals, Xsp$df, degree,1)
boundary_knots = c(0,1)

# function to run n integrals with the ibs method of the splines2 library
n_ibs_integrals <- function(n){
  for (i in 1:n){
    x = splines2::ibs(1, knots = internal_knots, degree = degree,
                  intercept = TRUE, Boundary.knots = boundary_knots) %*%
      coefs
    
  }
  return (x)
}

set.seed(43)
test_n <- c(1e0, 1e1, 1e2, .5e3, 1e3, .5e4)#, 1e4)
rationes2 <- numeric(length(test_n))

# run the microbenchmark
k = 1
for (i in test_n){
  a <- microbenchmark(
    "FdPot::compute_func_datum_integral" = 
      FdPot::compute_func_datum_integral(coefs, Xsp$argvals, Xsp$df, degree, 
                                         n_times = i),
    "splines2::ibs" = n_ibs_integrals(n = i)
  )
  a$time[a$expr == unique(a$expr)[1]]
  rationes2[k] = mean(a$time[a$expr != "splines2::ibs"]) / mean(a$time[a$expr == "splines2::ibs"])
  k = k + 1
}
# repeat the above code (I did it using the console)
x11()
par(mfrow=c(2,1))
plot(x=test_n, y=rationes1, col="navyblue", type="b", main="Efficiency of Rcpp integral, seed 42")
plot(x=test_n, y=rationes2, col="red", type="b", main="Efficiency of Rcpp integral, seed 43")

# test 2: use of openmp (to be run after compiling differently)
#====
# load data
# call microbenchmark NB CANNOT RUN SINCE R aborts it whenever num_threads > 1
#in the Makevars
#microbenchmark(
 # "compute_dissim_and_feats" = FdPot::test_case_compute_dissim_and_feats(X_coeffs=Xsp$fd$coefs,
  #                                          X_argvals=Xsp$argvals, 
 #                                           X_basis_df=Xsp$df,
#                                            X_basis_degree=degree, 
#                                            n_feats = 4)
#  )
#}
# run this manually 15 times
a <- numeric(15)

start.time <- Sys.time()
FdPot::test_case_compute_dissim_and_feats(X_coeffs=Xsp$fd$coefs,
                                          X_argvals=Xsp$argvals, 
                                          X_basis_df=Xsp$df,
                                          X_basis_degree=degree, 
                                          n_feats = 4
)
end.time <- Sys.time()
a[10] <- end.time - start.time

   

