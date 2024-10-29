
# cross-fitting folds
fold1 = sort(sample((m+1):n, round((n-m)/2)))
fold2 = setdiff((m+1):n, fold1)

fold1.m = sort(sample(1:m, round(m/2)))
fold2.m = setdiff(1:m, fold1.m)


#------------------------------------------
#    beta_delta
#------------------------------------------

h = sd(data$C[1:m])/(sum(data$delta[1:m]))^0.25
KC1 = dnorm(as.matrix(dist(data$C[fold1.m]/h,diag=T,upper=T)))/h
KC2 = dnorm(as.matrix(dist(data$C[fold2.m]/h,diag=T,upper=T)))/h

betadelta1 = init.beta(data$delta[fold1.m],data$Z[fold1.m,],KC1,
                       link = link, dlink = dlink)
betadelta2 = init.beta(data$delta[fold2.m],data$Z[fold2.m,],KC2,
                       link = link, dlink = dlink)

#------------------------------------------
#    beta_k
#------------------------------------------

m1 = length(fold1.m)
m2 = length(fold2.m)

betak1 <- foreach(iperturb = 1:Nperturb, .combine = cbind) %dopar%
  {
    V = rbeta(m1,0.5,1.5)*4
    out = init.beta.perturb(data$delta[fold1.m],data$Z[fold1.m,],KC1,
                                          V, 
                                          init = betadelta1,
                                          link = link, dlink = dlink)
    return(out)
  }

betak2 <- foreach(iperturb = 1:Nperturb, .combine = cbind) %dopar%
  {
    V = rbeta(m2,0.5,1.5)*4
    out = init.beta.perturb(data$delta[fold2.m],data$Z[fold2.m,],KC2,
                            V, 
                            init = betadelta2,
                            link = link, dlink = dlink)
    return(out)
  }


#------------------------------------------
#    Sk
#------------------------------------------

beta1.std = betadelta1/sqrt(sum(betadelta1^2))
beta2.std = betadelta2/sqrt(sum(betadelta2^2))
lp1 = drop(data$Z %*% beta1.std)
lp2 = drop(data$Z %*% beta2.std)
sdlp = sd(c(lp1,lp2))

# for (isetup in 1:Nsetup)
# {

h1 = sdlp/(sum(dataS$delta1))^0.3
h2 = sdlp/(sum(dataS$delta2))^0.3

Sk1 = rep(0,p*2)
Sk1[1:p] = Sk_sym(lp1[fold1],data$Z[fold1,],
                  dataS$X1[fold1],dataS$delta1[fold1],
                  data$C[fold1],dnorm,h1)
Sk1[p+1:p] = Sk_sym(lp1[fold1],data$Z[fold1,],
                    dataS$X2[fold1],dataS$delta2[fold1],
                    data$C[fold1],dnorm,h2)

Sk2 = rep(0,p*2)
Sk2[1:p] = Sk_sym(lp2[fold2],data$Z[fold2,],
                  dataS$X1[fold2],dataS$delta1[fold2],
                  data$C[fold2],dnorm,h1)
Sk2[p+1:p] = Sk_sym(lp2[fold2],data$Z[fold2,],
                    dataS$X2[fold2],dataS$delta2[fold2],
                    data$C[fold2],dnorm,h2)

#------------------------------------------
#    Perturbation 1: for Wk
#------------------------------------------


n1 = length(fold1)
n2 = length(fold2)

Skb1 <- foreach(iperturb = 1:Nperturb, .combine = rbind) %dopar%
  {
    V = matrix(rbeta(n1,0.5,1.5)*4)
   
    betak1.std = betak1[,iperturb]/sqrt(sum(betak1[,iperturb]^2))
    lp1 = (data$Z[fold1,] %*% betak1.std)
    out1 = drop(Sk_sym_perturb(lp1,data$Z[fold1,],
                               dataS$X1[fold1],dataS$delta1[fold1],
                               data$C[fold1],dnorm,h1,
                               V))
    out2 = drop(Sk_sym_perturb(lp1,data$Z[fold1,],
                               dataS$X2[fold1],dataS$delta2[fold1],
                               data$C[fold1],dnorm,h2,
                               V))
    return(c(out1,out2))
  }

Skb2 <- foreach(iperturb = 1:Nperturb, .combine = rbind) %dopar%
  {
    V = matrix(rbeta(n2,0.5,1.5)*4)
    
    betak2.std = betak2[,iperturb]/sqrt(sum(betak2[,iperturb]^2))
    lp2 = (data$Z[fold2,] %*% betak2.std)
    out1 = drop(Sk_sym_perturb(lp2,data$Z[fold2,],
                               dataS$X2[fold2],dataS$delta2[fold2],
                               data$C[fold2],dnorm,h2,
                               V))
    out2 = drop(Sk_sym_perturb(lp2,data$Z[fold2,],
                               dataS$X2[fold2],dataS$delta2[fold2],
                               data$C[fold2],dnorm,h2,
                               V))
    return(c(out1,out2))
  }

#------------------------------------------
#    Estimating W matrix
#------------------------------------------

W1.hat.nopen = W_hat_droplast_nopen(betak1, Skb1, npc)
W1.hat.ave = W_hat_droplast_ave(betak1, Skb1, betadelta1, npc)
W1.hat.ave.nopen = W_hat_droplast_ave_nopen(betak1, Skb1, betadelta1, npc)
W1.old = W_hat_perturb(betak1, Skb1)

W2.hat.nopen = W_hat_droplast_nopen(betak2, Skb2, npc)
W2.hat.ave = W_hat_droplast_ave(betak2, Skb2, betadelta2, npc)
W2.hat.ave.nopen = W_hat_droplast_ave_nopen(betak2, Skb2, betadelta2, npc)
W2.old = W_hat_perturb(betak2, Skb2)

# #------------------------------------------
# # CV for adaptive soft-thresholding SSL1
# #------------------------------------------
# 
# W1.hat = W1.hat.nopen
# W2.hat = W2.hat.nopen
# 
# betaSSL = (betadelta1 + betadelta2 
#            - drop(W1.hat %*% Sk2)-drop(W2.hat %*% Sk1))/2
# betaSSLk = (betak1 + betak2 
#             - W1.hat %*% t(Skb2) - W2.hat %*% t(Skb1))/2
# 
# betaSSL1 = betadelta1 - drop(W2.hat %*% Sk1)
# betaSSL2 = betadelta1 - drop(W2.hat %*% Sk1)
# 
# betaSSL.nopen = betaSSL
# betaSSLk.nopen = betaSSLk
# 
# #------------------------------------------
# # CV for adaptive soft-thresholding SSL2
# #------------------------------------------
# 
# W1.hat = W1.hat.ave
# W2.hat = W2.hat.ave
# 
# betaSSL = (betadelta1 + betadelta2 
#            - drop(W1.hat %*% Sk2)-drop(W2.hat %*% Sk1))/2
# betaSSLk = (betak1 + betak2 
#             - W1.hat %*% t(Skb2) - W2.hat %*% t(Skb1))/2
# 
# betaSSL.ave = betaSSL
# betaSSLk.ave = betaSSLk
# 
# #------------------------------------------
# # CV for adaptive soft-thresholding SSL3
# #------------------------------------------
# 
# W1.hat = W1.hat.ave.nopen
# W2.hat = W2.hat.ave.nopen
# 
# betaSSL = (betadelta1 + betadelta2 
#            - drop(W1.hat %*% Sk2)-drop(W2.hat %*% Sk1))/2
# betaSSLk = (betak1 + betak2 
#             - W1.hat %*% t(Skb2) - W2.hat %*% t(Skb1))/2
# 
# betaSSL.ave.nopen = betaSSL
# betaSSLk.ave.nopen = betaSSLk