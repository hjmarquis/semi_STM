dSk_asym = function(lp, Z, Xk, Dk, Ct, dK, h)
{
  n = nrow(Z)
  
  # ECDF of G
  X.order = order(Xk)
  C.sort = sort(Ct)
  GX = rep(0,n)
  Ctail = 0
  wZZ = matrix(0,n,n)
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
 
    Kbz.GXi2 = dK((lp[i]-lp[js])/h)/(h^2*GXi2)
    wZZ[i,i] = sum(Kbz.GXi2)
    diag(wZZ)[js] = diag(wZZ)[js] + Kbz.GXi2
    wZZ[i,js] = wZZ[js,i] = -Kbz.GXi2
  }
  
  return(drop(t(Z)%*% wZZ %*%Z)/w)
}