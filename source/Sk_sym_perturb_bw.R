Sk_sym_perturb_bw = function(lp, Z, Xk, Dk, Ct, K, h, V, 
                          max.bw = sqrt(-log(1e-8*sqrt(2*pi))))
{
  n = nrow(Z)
  n.perturb = ncol(V)
  lp = lp/h
  
  # Neighborhood info
  V.order = V.lw = matrix(1,n,n.perturb)
  V.up = matrix(n,n,n.perturb)
  for (i.perturb in 1:n.perturb) 
  {
    tmp.lp = lp[,i.perturb]
    tmp.order = order(tmp.lp)
    V.order[,i.perturb] = tmp.order
    tmp.lw = 1
    i.lw = tmp.order[tmp.lw]
    for (tmp.up in 2:n)
    {
      i.up = tmp.order[tmp.up]
      if(tmp.lp[i.up] - tmp.lp[i.lw] > max.bw)
      {
        V.up[i.lw,i.perturb] = tmp.up-1
        for(tmp.lw in (tmp.lw+1):tmp.up)
        {
          i.lw = tmp.order[tmp.lw]
          if(tmp.lp[i.up] - tmp.lp[i.lw] <= max.bw)
            break
          V.up[i.lw,i.perturb] = tmp.up-1
        }
      }
      V.lw[i.up,i.perturb] = tmp.lw
    }
  }
  
  # ECDF of G
  X.order = order(Xk)
  C.order = order(Ct)
  GX = rep(0,n)
  Ctail.count = 0
  Ctail.sum = rep(0,n.perturb)
  wZ = matrix(0,n.perturb,n)
  w = rep(0,n.perturb)
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
      Ctail.sum = Ctail.sum + V[C.order[n-Ctail.count],]
      Ctail.count = Ctail.count + 1
    }
    
    if( (next.X == n) | (Dk[i]== 0))
      next
    GXi2 = Ctail.sum^2    
    js = X.order[(next.X+1):n]
    njs = n-next.X
    
    for (i.perturb in 1:n.perturb)
    {
      w[i.perturb] = (w[i.perturb] + 
        sum(V[i,i.perturb]*V[js,i.perturb]/GXi2[i.perturb]))
      ijs.neibor = V.order[
        V.lw[i,i.perturb]:V.up[i,i.perturb]
        ,i.perturb]
      ijs = intersect(ijs.neibor,js)
      if(length(ijs)==0)
        next
      Vij.GXi2 = V[i,i.perturb]*V[ijs,i.perturb]/GXi2[i.perturb]
      
      Kbz.GXi2 = K(lp[i,i.perturb]-lp[ijs,i.perturb])*Vij.GXi2/h 
      wZ[i.perturb,i] = sum(Kbz.GXi2)
      wZ[i.perturb,ijs] = wZ[i.perturb,ijs] - Kbz.GXi2
    }
  }
  
  return(drop(wZ%*%Z)/w)
}