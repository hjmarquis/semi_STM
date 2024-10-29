ddSk_sym = function(lp, k, Z, Xk, Dk, Ct, ddK, h)
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
    GXi2 = (Ctail/n)^2 
    w = w + (n-next.X)/GXi2
   
    js = X.order[(next.X+1):n]
    Kbz.GXi2 = (Z[i,k]-Z[js,k])*ddK((lp[i]-lp[js])/h)/(h^3*GXi2)
    wZZ[i,i] = sum(Kbz.GXi2)
    wZZ[diag.pos[js]] = wZZ[diag.pos[js]] + Kbz.GXi2
    wZZ[i,js] = wZZ[js,i] = -Kbz.GXi2
  }
  
  return((t(Z)%*% wZZ %*%Z)/w)
}