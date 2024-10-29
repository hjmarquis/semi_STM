cv.concordance = function(lp,Xk, Dk, Ct)
{
  n = nrow(lp)
  nlam = ncol(lp)
  
  # ECDF of G
  X.order = order(Xk)
  C.sort = sort(Ct)
  Ctail = 0
  wc = rep(0,nlam)
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
    wc = wc + apply(lp[i,] > t(lp[js,]),1,sum)/GXi2
  }
  
  return(cbind(wc,w))
}