Sk_asym = function(lp, Z, Xk, Dk, Ct, K, h)
{
  n = nrow(Z)
  
  # ECDF of G
  X.order = order(Xk)
  C.sort = sort(Ct)
  GX = rep(0,n)
  Ctail = 0
  wZ = rep(0,n)
  w = 0
  next.X = n
  for(i.X in (n-1):1)
  {
    i = X.order[i.X]
    if(Xk[i]<Xk[X.order[next.X]])
      next.X = i.X
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
    w = w + njs/GXi2
 
    Kbz = K((lp[i]-lp[js])/h)/h
    wZ[i] = sum(Kbz)/GXi2
    wZ[js] = wZ[js] - Kbz/GXi2
  }
  
  return(drop(wZ%*%Z)/w)
}