source("dlm_ar_functions.r")
require(dlm)

# Set data and graphics path
dpath = "/storage/sheinson_research/"
gpath = "../graphs/"

dlm.mle <- function(nt, nsims, mod, x, beta0 = 750, beta1 = 15, phi = 0.6, rho = 0.8, sigma2s = 1, sigma2b = 1, sigma2m = 1, print.iter=FALSE)
{
  # Create U matrix
  if(x == "conv")
  {
    load(paste(dpath,"fmri-design-id.rdata",sep=""))
    ind = which(sapply(fmri.design, function(a) length(a$t)) == nt)
    U = fmri.design[[ind]]$X
  } else if(x == "norm") {
    U = cbind(1,rnorm(nt))
  } else if(x == "none") {
    U = cbind(1,rep(0,nt))
  }
  
  # Create FF matrix
  if(mod == "M011")
  {
    if(x == "none") FF = U[,1,drop=FALSE] else FF = U[,2,drop=FALSE]
  } else if(mod == "M101") {
    FF = U[,1,drop=FALSE]
  } else {
    FF = U
  }
  
  # Create storage for MLEs
  theta.mle = matrix(NA, nr = nsims, nc = 7)
  error = numeric(nsims)
  converge = rep(NA, nsims)
  colnames(theta.mle) = c("beta0", "beta1", "phi", "rho", "sigma2s", "sigma2b","sigma2m")
  
  # Start loop over simulations
  for(i in 1:nsims)
  {
    # Simulate time series
    if(mod == "M011" | mod == "M101")
    {
      sim = dlm.sim(nt, c(beta0,beta1), phi, sigma2s, sigma2m, U, FF)
    } else {
      sim = dlm.sim(nt, c(beta0,beta1), c(phi,rho), c(sigma2s, sigma2b), sigma2m, U, FF)
    }
    
    # Calculate MLE
    if(x == "none") U.ls = U[,1,drop=FALSE] else U.ls = U
    fit0 = lm(sim$y ~ U.ls - 1)
    phi.init = as.numeric(acf(residuals(fit0),type='partial',plot=FALSE)$acf)[1]
    sigma2s.init = summary(fit0)$sigma^2
    sigma2m.init = summary(fit0)$sigma^2
    if(mod != "M011" & mod != "M101")
    {
      phi.init = c(phi.init,phi.init)
      sigma2s.init = c(sigma2s.init/2,sigma2s.init/2)
    }
    fit <- try(dlmMLE(sim$y, c(phi.init,log(sigma2s.init),log(sigma2m.init)), function(par) build.ar(par, U, FF), method="Nelder-Mead"), silent=TRUE)
    if(class(fit) != 'try-error')
    {
      converge[i] = fit$convergence
      s <- dlmSmooth(dlmFilter(sim$y, build.ar(fit$par, U, FF)))
      if(mod == "M011" | mod == "M101")
      {
        theta.mle[i,c(1,2,3,5,7)] <- c(s$s[nt+1,1:2],fit$par[1],exp(fit$par[2:3]))
      } else{
        theta.mle[i,] <- c(s$s[nt+1,1:2],fit$par[1:2],exp(fit$par[3:5]))
      }
    } else {
      error[i] = 1
    }
    if(print.iter) print(paste(nt, mod, nsims, x, i, error[i], converge[i]))
  }
  
  file=paste(dpath,paste(nt,nsims,mod,x,beta0,beta1,phi,rho,sigma2s,sigma2b,sigma2m,sep="-"),".rdata",sep="")
  print(file)
  out.mle = list(theta.mle=theta.mle,error=error,converge=converge)
  save(out.mle, file=file)
}

require(plyr)
require(doMC)
registerDoMC()
#data1 = expand.grid(nsims = 1000, mod=c("M011","M101"), beta0 = 750, beta1 = 15, phi = c(.1,.5,.9), rho=.5, sigma2s = c(1,5,10,15,20), sigma2b = 10, sigma2m = 10, nt = c(250,500,1000), x = "conv", stringsAsFactors=FALSE)
data1 = expand.grid(nsims = 1000, mod="M111", beta0 = 750, beta1 = 15, phi = c(.3,.6,.9), rho = c(.3,.6,.9), sigma2s = c(1,5,10,15,20), sigma2b = c(1,5,10,15,20), sigma2m = 10, nt = c(250,500,1000), x = "conv", stringsAsFactors=FALSE)
data2 = expand.grid(nsims = 1000, mod = c("M011","M101"), beta0 = 750, beta1 = 15, phi = 0.5, rho=.5, sigma2s = 5, sigma2b = 10, sigma2m = 10, nt = c(250,500,1000), x = "norm", stringsAsFactors=FALSE)
data3 = expand.grid(nsims = 1000, mod = "M111", beta0 = 750, beta1 = 15, phi = 0.6, rho=0.6, sigma2s = 10, sigma2b = 10, sigma2m = 10, nt = c(250,500,1000), x = "norm", stringsAsFactors=FALSE)
data4 = expand.grid(nsims = 1000, mod = "M101", beta0 = 750, beta1 = 15, phi = 0.5, rho=0.5, sigma2s = 5, sigma2b = 10, sigma2m = 10, nt = c(250,500,1000), x = "none", stringsAsFactors=FALSE)
mydata = rbind(data1,data2,data3,data4)
m_ply(mydata, dlm.mle, .parallel=TRUE)
