# Function to simulate nt + 1 univariate states (x_t) and nt univariate observations (y_t) from an AR(1) dlm
# y_t = U_t*beta + F_t*x_t + v_t, v_t iid N(0,sigma2m)
# x_t = phi*x_t-1 + w_t, w_t iid N(0,sigma2s)
# U is an nt by d matrix of regression covariates
# FF is an nt length vector of coefficients for the state
dlm.sim <- function(nt, beta, phi, sigma2s, sigma2m, U, FF, x0=0)
{
  y = numeric(nt)
  x = matrix(NA, nr = length(phi), nc= nt+1)
  x[,1] = x0
  v = rnorm(nt, 0, sqrt(sigma2m))
  w = matrix(NA, nr = length(phi), nc= nt)
  for(i in 1:length(phi)) w[i,] = rnorm(nt, 0, sqrt(sigma2s[i]))
  for(i in 1:nt)
  {
    x[,i+1] = phi*x[,i] + w[,i]
    y[i] = t(U[i,])%*%beta + FF[i,]%*%x[,i+1] + v[i]
  }
  return(list(y=y,x=x,true.params=list(theta=c(beta,phi,sigma2s,sigma2m),U=U,FF=FF)))
}

# Function to build a 'dlm' object for a regression model with AR(1) states
# No constraint on autocorrelation parameter phi
# Assume univariate data
# par = c(phi1, phi2, ..., phip, log(sigma2s1), log(sigma2s2), ..., log(sigma2sp), log(sigma2m))
# U is an nt by d matrix of regression covariates
# f is an nt by p matrix of coefficients for the state

build.ar <- function(par, U, f)
{
  require(dlm)
  d = dim(U)[2]
  p = dim(f)[2]
  nt = dim(U)[1]
  FF = matrix(rep(0,d+p),nr=1)
  JFF = matrix(1:(d+p),nr=1)
  X = cbind(U,f)
  V = exp(par[2*p + 1])
  GG = diag(d+p)
  for(i in 1:p) GG[d+i,d+i] = par[i]
  W = 0*diag(d+p)
  for(i in 1:p) W[d+i,d+i] = exp(par[p + i])
  m0 = rep(0,d+p)
  C0 = 1e6*diag(d+p); C0[(d+1):(d+p),(d+1):(d+p)] = 0
  mod = dlm(FF=FF,V=matrix(V),GG=GG,W=W,m0=m0,C0=C0,JFF=JFF,X=X)
  return(mod)
}
