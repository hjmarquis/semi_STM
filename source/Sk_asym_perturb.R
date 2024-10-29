Sk_asym_perturb = function(lp, Z, Xk, Dk, Ct, K, h, V)
{
  n = nrow(Z)
  n.perturb = ncol(V)
  
  # ECDF of G
  X.order = order(Xk)
  C.sort = sort(Ct)
  GX = rep(0,n)
  Ctail = 0
  wZ = matrix(0,n.perturb,n)
  w = rep(0,n.perturb)
  next.X = n
  for(i.X in (n-1):1)
  {
    i = X.order[i.X]
    if(Xk[i]<Xk[X.order[next.X]])
      next.X = i.X
    if(Dk[i]==0)
      next
    while (Ctail < n) 
    {
      if(C.sort[n-Ctail] < Xk[i])
        break
      Ctail = Ctail + 1
    }
    GXi2 = (Ctail/n)^2  
    js = intersect(X.order[next.X:n],(i+1):n)
    njs = length(js)
    
    if( (next.X == n) | (Dk[i]== 0) |
        njs==0)
      next  
    
    Vij = V[i,]*matrix(t(V[js,]),n.perturb)
    w = w + apply(Vij,1,sum)/GXi2
    
    Kbz = K((lp[i,]-
               matrix(t(lp[js,]),n.perturb))/h)/h * Vij
    wZ[,i] = apply(Kbz,1,sum)/GXi2
    wZ[,js] = wZ[,js] - Kbz/GXi2
  }
  
  return(drop(wZ%*%Z)/w)
}