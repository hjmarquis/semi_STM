dSk_sym_perturb = function(lp, Z, Xk, Dk, Ct, dK, h, V)
{
  n = nrow(Z)
  
  # ECDF of G
  diag.pos = n*(1:n - 1) + 1:n
  X.order = order(Xk)
  C.order = order(Ct)
  GX = rep(0,n)
  Ctail.count = 0
  Ctail.sum = 0
  wZZ = matrix(0,n,n)
  w = 0
  next.X = n
  for(i.X in (n-1):1)
  {
    i = X.order[i.X]
    if(Xk[i]<Xk[X.order[next.X]])
      next.X = i.X
    while (Ctail.count < n) 
    {
      if(Ct[C.order[n-Ctail.count]] < Xk[i])
        break
      Ctail.sum = Ctail.sum + V[C.order[n-Ctail.count]]
      Ctail.count = Ctail.count + 1
    }
    if( (next.X == n) | (Dk[i]== 0))
      next
    GXi2 = Ctail.sum^2
   
    js = X.order[(next.X+1):n]
    
    Vij.GXi2 = V[i]*t(V[js])/GXi2
    w = w + sum(Vij.GXi2)
    Kbz.GXi2 = dK((lp[i]-lp[js])/h)/(h^2)*Vij.GXi2
    
    wZZ[i,i] = sum(Kbz.GXi2)
    wZZ[diag.pos[js]] = wZZ[diag.pos[js]] + Kbz.GXi2
    wZZ[i,js] = wZZ[js,i] = -Kbz.GXi2
  }
  
  return((t(Z)%*% wZZ %*%Z)/w)
}