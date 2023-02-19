library(FdPot)
# example from fda
data_W <- CanadianWeather$dailyAv[,,1]
m <- 5           # spline order 
degree <- m-1    # spline degree 
nbasis = 20
basis <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis, norder=m)
time = seq(0, 1, length.out = 365)
Xsp <- smooth.basis(argvals=time, y=data_W, fdParobj=basis)
Xsp$fd$basis$basisvalues

# now use the function
b <- rcpp_basis(time, df = Xsp$df, degree = degree)
b.02 <- rcpp_basis_constr_due(c(0.5, .06), x_orig = time,  df = Xsp$df, degree = degree)
dim(b.02$basismat)
eval.fd(time, Xsp$fd)
dim(b$basismat)
Xsp$fd$basis$nbasis

b$basismat %*% Xsp$fd$coefs
b
# evaluate test 1
all.equal(b$basismat %*% Xsp$fd$coefs, eval.fd(time, Xsp$fd))
all.equal(b.02$basismat %*% Xsp$fd$coefs, eval.fd(time, Xsp$fd))



# test 2: eval on single point
set.seed(20)
idx <- which(time == sample(time, size = 1))
b2 <- rcpp_basis(time[idx:(idx+1)], df = Xsp$df, degree = degree)
# TODO error here.
b3  <- rcpp_basis(time[c(1, idx, 365)], df = Xsp$df, degree = degree)
# b2. NOT A MATCH
all.equal(b2$basismat %*% Xsp$fd$coefs, eval.fd(matrix(time[idx:(idx+1)], ncol=1), Xsp$fd))
# b restricted to point of interest
all.equal(b$basismat[idx,] %*% Xsp$fd$coefs, eval.fd(matrix(time[idx], ncol=1), Xsp$fd))
# b3 (pick elt in the middle)
all.equal(b3$basismat[2,] %*% Xsp$fd$coefs,  eval.fd(matrix(time[idx], ncol=1), Xsp$fd))
# Conclusion: it does not work because the constructor uses the input vector x to construct the
# basis.

# use b.03
all.equal(b$basismat[idx,] %*% Xsp$fd$coefs, eval.fd(matrix(time[idx], ncol=1), Xsp$fd))
new_x <- c(.533, .534)
all.equal(b.02$basismat %*% Xsp$fd$coefs, eval.fd(matrix(new_x,ncol=1), Xsp$fd))
      