dSk_sym_fix = function(lp, Z, Xk, Dk, Ct, dK, h)
{
  n = nrow(Z)
  
  # ECDF of G
  diag.pos = n*(1:n - 1) + 1:n
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
    if( (next.X == n) | (Dk[i]== 0))
      next
    GXi = (Ctail/n)
    w = w + (n-next.X)/GXi
   
    js = X.order[(next.X+1):n]
    Kbz.GXi = dK((lp[i]-lp[js])/h)/(h^2*GXi)
    wZZ[i,i] = sum(Kbz.GXi)
    wZZ[diag.pos[js]] = wZZ[diag.pos[js]] + Kbz.GXi
    wZZ[i,js] = wZZ[js,i] = -Kbz.GXi
  }
  
  return(drop(t(Z)%*% wZZ %*%Z)/w)
}