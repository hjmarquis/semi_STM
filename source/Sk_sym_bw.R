Sk_sym_bw = function(lp, Z, Xk, Dk, Ct, K, h,
                     max.bw = sqrt(-log(1e-8*sqrt(2*pi))))
{
  n = nrow(Z)
  n.perturb = 1
  lp = lp/h
  
  # Neighborhood info
  V.order = V.lw = rep(1,n)
  V.up = rep(n,n)

  tmp.lp = lp
  V.order = tmp.order = order(lp)
  tmp.lw = 1
  i.lw = tmp.order[tmp.lw]
  for (tmp.up in 2:n)
  {
    i.up = tmp.order[tmp.up]
    if(tmp.lp[i.up] - tmp.lp[i.lw] > max.bw)
    {
      V.up[i.lw] = tmp.up-1
      for(tmp.lw in (tmp.lw+1):tmp.up)
      {
        i.lw = tmp.order[tmp.lw]
        if(tmp.lp[i.up] - tmp.lp[i.lw] <= max.bw)
          break
        V.up[i.lw] = tmp.up-1
      }
    }
    V.lw[i.up] = tmp.lw
  }
  
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
    
    if( (next.X == n) | (Dk[i]== 0))
      next    
    js = X.order[(next.X+1):n]
    w = w + (n-next.X)/GXi2
    

    js.neibor = V.order[
      V.lw[i]:V.up[i]
      ]
    ijs = intersect(js.neibor,js)
    if(length(ijs)==0)
      next
    
    Kbz.GXi2 = K(lp[i]-lp[ijs])/(GXi2*h) 
    wZ[i] = sum(Kbz.GXi2)
    wZ[ijs] = wZ[ijs] - Kbz.GXi2

  }
  
  return(drop(wZ%*%Z)/w)
}