Sk_sym_perturb_bw_fix = function(lp, Z, Xk, Dk, Ct, K, h, V, 
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
    V.order[,i.perturb] = tmp.order = order(lp[,i.perturb])
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
  C.sort = sort(Ct)
  Ctail = 0
  wZ = matrix(0,n.perturb,n)
  w = rep(0,n.perturb)
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
    GXi = (Ctail/n)
    
    if( (next.X == n) | (Dk[i]== 0))
      next    
    js = X.order[(next.X+1):n]
    njs = n-next.X
    
    for (i.perturb in 1:n.perturb)
    {
      w[i.perturb] = (w[i.perturb] + 
        sum(V[i,i.perturb]*V[js,i.perturb]/GXi))
      js.neibor = V.order[
        V.lw[i,i.perturb]:V.up[i,i.perturb]
        ,i.perturb]
      js = intersect(js.neibor,js)
      if(length(js)==0)
        next
      Vij.GXi = V[i,i.perturb]*V[js,i.perturb]/GXi
      
      Kbz.GXi = K(lp[i,i.perturb]-lp[js,i.perturb])*Vij.GXi/h 
      wZ[i.perturb,i] = sum(Kbz.GXi)
      wZ[i.perturb,js] = wZ[i.perturb,js] - Kbz.GXi
    }
  }
  
  return(drop(wZ%*%Z)/w)
}