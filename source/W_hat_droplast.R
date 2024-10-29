require(glmnet)
W_hat_droplast = function(betak,Skb,npc = nrow(betak)-1)
{
  p = nrow(betak)
  B = ncol(betak)
  Kp = ncol(Skb)
  K = Kp/p
  Proj = diag(1,Kp)[,setdiff(1:Kp,(1:K-1)*p+1)[1:npc]]
  PSkb =  Skb %*% Proj
  
  W.hat = matrix(0,p,npc)
  
  for(j in 1:p)
  {
    tmp.ridge = cv.glmnet(PSkb,betak[j,],alpha = 0)
    W.hat[j,]=coef(tmp.ridge,s="lambda.min")[-1]
  }
  
  return(W.hat %*% t(Proj))
}