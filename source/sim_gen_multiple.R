sim.gen.multiple = function(n,p,beta0,sigmaZ,mu,sigma,m, Nperturb=0,
                            inv.link = logit, inv.h = exp)
{
  Z = mvrnorm(n,rep(0,p),sigmaZ^2*(0.2+0.8*diag(1,p)))
  u = runif(n)
  T = inv.h((inv.link(u) - Z%*%beta0)/3)*4 #h^-1 (g^-1(u) - beta'Z)
  C = runif(n,0,12)
  delta = (T<=C)
  if(m < n)
    delta[(m+1):n] = NA
  
  nsetup = length(mu)
  error = matrix(rnorm(n*4),n,4)
  
  Di = rbinom(n,1,0.5)  
  Di2 = rbinom(n,1,0.5)
  
  dataS = vector("list",nsetup)
  for(isetup in 1:nsetup)
  {
    epsilon = (Di*(error[,1]*sigma[[isetup]][1]+mu[[isetup]][1]) +
                (1-Di)*(error[,2]*sigma[[isetup]][2]+mu[[isetup]][2]))
    Tstar = exp(epsilon)*T
    X1 = pmin(Tstar,C)
    delta1 = (Tstar<=C)
  
    epsilon2 = (Di2*(error[,3]*sigma[[isetup]][3]+mu[[isetup]][3]) +
                  (1-Di2)*(error[,4]*sigma[[isetup]][4]+mu[[isetup]][4]))
    Tstar2 = exp(epsilon2)*T
    X2 = pmin(Tstar2,C)
    delta2 = (Tstar2<=C)
    
    dataS[[isetup]] = data.frame(X1=X1,delta1=delta1,
                                 X2=X2,delta2=delta2)
  }
  
  data = data.frame(delta = delta, C=C)
  data$Z = Z
  if(Nperturb > 0)
    data$V = matrix(rbeta(n*Nperturb,0.5,1.5)*4,n)
  
  return(list(data=data,dataS=dataS))
}
