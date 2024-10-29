concordance12 = function(X1,D1,X2,D2, G1, G2)
{
  n = length(X1)
  
  # ECDF of G
  X.order = order(X1)
  wc = 0
  w = 0
  next.X = n
  for(i.X in (n-1):1)
  {
    i = X.order[i.X]
    if(X1[i]<X1[X.order[next.X]])
      next.X = i.X
    
    if( (next.X == n) | (D1[i]== 0))
      next
    js = X.order[(next.X+1):n]
    
    js0 = js[((X2[js]< X2[i])& (D2[js] == 1))]
    js1 = js[((X2[js]> X2[i])& (D2[i] == 1))]
    
    w1 = length(js1)/max(G1[i],G2[i])
    w0 = ifelse(length(js0)==0, 0, 
                sum(1/pmax(G2[js0],G1[i])))
    
    w = w + w1 + w0
    
    wc = wc + w1
  }
  
  return(wc/w)
}