sim.gen = function(n,p,beta0,sigmaZ,mu,sigma,m, Nperturb=0)
{
  Z = mvrnorm(n,rep(0,p),sigmaZ^2*(0.2+0.8*diag(1,p)))
  u = runif(n)
  T = exp((log(u/(1-u)) - Z%*%beta0)/3)*4 #h^-1 (g^-1(u) - beta'Z)
  C = runif(n,0,12)
  delta = (T<=C)
  if(m < n)
    delta[(m+1):n] = NA
  
  Di = rbinom(n,1,0.5)
  epsilon = Di*rnorm(n,mu[1],sigma[1]) + (1-Di)*rnorm(n,mu[2],sigma[2])
  Tstar = exp(epsilon)*T
  X1 = pmin(Tstar,C)
  delta1 = (Tstar<=C)
  
  Di2 = rbinom(n,1,0.5)
  epsilon2 = Di2*rnorm(n,mu[3],sigma[3]) + (1-Di2)*rnorm(n,mu[4],sigma[4])
  Tstar2 = exp(epsilon2)*T
  X2 = pmin(Tstar2,C)
  delta2 = (Tstar2<=C)
  
  out = data.frame(delta = delta, C=C, 
                   X1=X1,delta1=delta1,
                   X2=X2,delta2=delta2)
  out$Z = Z
  if(Nperturb > 0)
    out$V = matrix(rbeta(n*Nperturb,0.5,1.5)*4,n)
  
  return(out)
}