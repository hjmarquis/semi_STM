require(glmnet)
W_hat_perturb_no_intercept = function(betak,Skb)
{
  p = nrow(betak)
  Kp = ncol(Skb)
  W.hat = matrix(0,p,Kp)
  
  for(j in 1:p)
  {
    tmp.ridge = cv.glmnet(Skb,betak[j,],intercept=FALSE,alpha = 0)
    W.hat[j,]=coef(tmp.ridge,s="lambda.min")[-1]
  }
  
  return(W.hat)
}