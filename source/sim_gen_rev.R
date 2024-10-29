sim.gen.rev = function(n,beta0,sigmaZ,m,
                       misT, Nperturb=0,
                       inv.link = logit, 
                       inv.h = function(x) exp(x/3)*4)
{
  p = length(beta0)
  Z = mvrnorm(n,rep(0,p),sigmaZ^2*(0.2+0.8*diag(1,p)))
  u = runif(n)
  Tevent = inv.h((inv.link(u) - Z%*%beta0)) #h^-1 (g^-1(u) - beta'Z)
  C = runif(n,0,12)
  delta = (Tevent<=C)
  if(m < n)
    delta[(m+1):n] = NA
  
  Tstar = misT(Tevent)

  X1 = pmin(Tstar[,1],C)
  delta1 = (Tstar[,1]<=C)

  X2 = pmin(Tstar[,2],C)
  delta2 = (Tstar[,2]<=C)
  
  dataS = data.frame(X1=X1,delta1=delta1,
                     X2=X2,delta2=delta2)
  
  data = data.frame(delta = delta, C=C,
                    Tevent = Tevent, Tstar = Tstar)
  data$Z = Z
  if(Nperturb > 0)
    data$V = matrix(rbeta(n*Nperturb,0.5,1.5)*4,n)
  
  return(list(data=data,dataS=dataS))
}

noise.loglinear.S = function(t)
{
  mu = c(0,0.5,-0.25,0)
  sigma = c(0.5,0.15,0.35,0.45)
  
  n = length(t)
  
  Di = rbinom(n,1,0.5)  
  Di2 = rbinom(n,1,0.5)
  
  error = matrix(rnorm(n*4),n,4)
  
  epsilon = (Di*(error[,1]*sigma[1]+mu[1]) +
               (1-Di)*(error[,2]*sigma[2]+mu[2]))
  Tstar = exp(epsilon)*t
  
  epsilon2 = (Di2*(error[,3]*sigma[3]+mu[3]) +
                (1-Di2)*(error[,4]*sigma[4]+mu[4]))
  Tstar2 = exp(epsilon2)*t
  
  return(cbind(Tstar,Tstar2))
}

noise.loglinear.L = function(t)
{
  mu = c(1,-0.5,0,1.5)
  sigma = c(1.5,0.5,1,0.5)
  
  n = length(t)
  
  Di = rbinom(n,1,0.5)  
  Di2 = rbinom(n,1,0.5)
  
  error = matrix(rnorm(n*4),n,4)
  
  epsilon = (Di*(error[,1]*sigma[1]+mu[1]) +
               (1-Di)*(error[,2]*sigma[2]+mu[2]))
  Tstar = exp(epsilon)*t
  
  epsilon2 = (Di2*(error[,3]*sigma[3]+mu[3]) +
                (1-Di2)*(error[,4]*sigma[4]+mu[4]))
  Tstar2 = exp(epsilon2)*t
  
  return(cbind(Tstar,Tstar2))
}

noise.exp.S = function(t)
{
  n = length(t)
  r1 = 0.9
  r2 = 0.05
  
  Tstar = pmax(t+ r2*(rexp(n)-0.5),0)
  Tstar2 = pmax(t+ r2*(rchisq(n,df=1)-0.5),0)
  
  return(cbind(Tstar,Tstar2))
}

noise.exp.L = function(t)
{
  n = length(t)
  r1 = 0.45
  r2 = 0.5
  
  Tstar = pmax(t+ r2*(rexp(n)-0.5),0)
  Tstar2 = pmax(t+ r2*(rchisq(n,df=1)-0.5),0)
  
  return(cbind(Tstar,Tstar2))
}
