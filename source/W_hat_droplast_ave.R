require(glmnet)
W_hat_droplast_ave = function(betak,Skb, beta.thres,npc = nrow(betak)-1)
{
  p = nrow(betak)
  B = ncol(betak)
  Kp = ncol(Skb)
  K = Kp/p
  
  k.nz = which(beta.thres!=0)
  n.ave = max(1,length(k.nz))
  
  if(npc == nrow(betak)-1)
  {
    proj.list = outer(1:p, p*rep(1:K-1,each=n.ave), '+')
    proj.list[k.nz+p*(1:(n.ave*K)-1)] = 0
  }else
  {
    proj.list = matrix(1:Kp,Kp, n.ave)
    proj.list[rep(k.nz,each=K)+p*(1:(n.ave*K)-1)] = 0
  }
  
  W.hat.ave = 0
  for(k in 1:ncol(proj.list))
  {
    Proj = diag(1,Kp)[,proj.list[,k]]
    PSkb =  Skb %*% Proj
    
    W.hat = matrix(0,p,npc)
    
    for(j in 1:p)
    {
      tmp.ridge = cv.glmnet(PSkb,betak[j,],alpha = 0)
      W.hat[j,]=coef(tmp.ridge,s="lambda.min")[-1]
    }
    
    W.hat.ave = W.hat.ave + W.hat %*% t(Proj)/ncol(proj.list)
  }
  
  return(W.hat.ave)
}