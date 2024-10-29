baseSTM = function(delta,Z, KC, 
                  V, init = rep(0,ncol(Z)), tol=1e-7,
                  maxit = 100,min.factor = 0.75,
                  ls.factor = 0.75,max.move = 1,
                  link = expit, dlink = dexpit)
{
  n = nrow(Z)
  KC = t(V*KC)
  KCd = drop(KC%*%delta)
  hC = rep(0,n)
  oldscore = NULL
  max.dbeta = max.move
  
  # print(k)
  # print(init)
  # print(oldscore)
  lp = drop(Z%*%init)
  hC.flag = rep(TRUE,n)
  gij = wZbar = matrix(0,n,n)
  hHess = rep(0,n)
  for(kk in 1:maxit)
  {
    lij = outer(hC[hC.flag],lp,"+")
    gij[hC.flag,] = link(lij)
    tmp = KC[hC.flag,]*gij[hC.flag,]
    wZbar[hC.flag,] = KC[hC.flag,]*dlink(lij)
    if(sum(hC.flag)>=2)
    {
      hscore = apply(tmp,1,sum)-KCd[hC.flag]
      hHess[hC.flag] = apply(wZbar[hC.flag,],1,sum)
    }else
    {
      hscore = sum(tmp)-KCd[hC.flag]
      hHess[hC.flag] = sum(wZbar[hC.flag,])
    }
    
    dhC = hscore/hHess[hC.flag]
    dhC = sign(dhC)*pmin(abs(dhC),max.move)
    kk.flag = abs(hscore) > tol
    if(!any(kk.flag))
      break
    hC[hC.flag][kk.flag] = hC[hC.flag][kk.flag] - dhC[kk.flag]
    hC.flag[hC.flag] = kk.flag
  }
  if(kk >= maxit)
    stop("Numerical error when computing h0(Ci)")
  return(hC)
}