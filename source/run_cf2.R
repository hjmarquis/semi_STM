# cross-fitting folds
fold1 = sort(sample((m+1):n, round((n-m)/2)))
fold2 = setdiff((m+1):n, fold1)

# fold1.m = sort(sample(1:m, round(m/2)))
# fold2.m = setdiff(1:m, fold1.m)

#------------------------------------------
#    beta_delta
#------------------------------------------

h = sd(data$C[1:m])/(sum(data$delta[1:m]))^0.25
KC = dnorm(as.matrix(dist(data$C[1:m]/h,diag=T,upper=T)))/h

betadelta = init.beta(data$delta[1:m],data$Z[1:m,],KC,
                       link = link, dlink = dlink)

#------------------------------------------
#    Sk
#------------------------------------------

beta.std = betadelta/sqrt(sum(betadelta^2))
lp = drop(data$Z %*% beta.std)
sdlp = sd(lp)

# for (isetup in 1:Nsetup)
# {

h1 = sdlp/(sum(dataS$delta1))^0.3
h2 = sdlp/(sum(dataS$delta2))^0.3
Sk1 = rep(0,p*2)
Sk1[1:p] = Sk_sym(lp[fold1],data$Z[fold1,],
                  dataS$X1[fold1],dataS$delta1[fold1],
                  data$C[fold1],dnorm,h1)
Sk1[p+1:p] = Sk_sym(lp[fold1],data$Z[fold1,],
                    dataS$X2[fold1],dataS$delta2[fold1],
                    data$C[fold1],dnorm,h2)

Sk2 = rep(0,p*2)
Sk2[1:p] = Sk_sym(lp[fold2],data$Z[fold2,],
                  dataS$X1[fold2],dataS$delta1[fold2],
                  data$C[fold2],dnorm,h1)


#------------------------------------------
#    Perturbation 1: for Wk
#------------------------------------------

betak1 =  matrix(0,p, Nperturb)    
betak2 =  matrix(0,p, Nperturb)
Skb1 = matrix(0, Nperturb,p*2)
Skb2 = matrix(0, Nperturb,p*2)     
for (iperturb in 1:Nperturb) 
{
  V = rbeta(n,0.5,1.5)*4
  betak[,iperturb] = init.beta.perturb(data$delta[1:m],data$Z[1:m,],KC1,
                                        V[1:m], 
                                        init = betadelta,
                                        link = link, dlink = dlink)
  
  betak.std = betak[,iperturb]/sqrt(sum(betak[,iperturb]^2))
  lp = (data$Z %*% betak.std)
  
  Skb1[iperturb,1:p] = drop(Sk_sym_perturb(lp[fold1,,drop=F],data$Z[fold1,],
                                           dataS$X1[fold1],dataS$delta1[fold1],
                                           data$C[fold1],dnorm,h1,
                                           matrix(V[fold1])))
  Skb1[iperturb,p+1:p] = drop(Sk_sym_perturb(lp[fold1,,drop=F],data$Z[fold1,],
                                             dataS$X2[fold1],dataS$delta2[fold1],
                                             data$C[fold1],dnorm,h2,
                                             matrix(V[fold1])))  
  
  Skb2[iperturb,1:p] = drop(Sk_sym_perturb(lp[fold2,,drop=F],data$Z[fold2,],
                                           dataS$X1[fold2],dataS$delta1[fold2],
                                           data$C[fold2],dnorm,h1,
                                           matrix(V[fold2])))
  Skb2[iperturb,p+1:p] = drop(Sk_sym_perturb(lp[fold2,,drop=F],data$Z[fold2,],
                                             dataS$X2[fold2],dataS$delta2[fold2],
                                             data$C[fold2],dnorm,h2,
                                             matrix(V[fold2])))  
}

#------------------------------------------
#    Estimating W matrix
#------------------------------------------

W1.hat.nopen = W_hat_droplast_nopen(betak, Skb1, npc)
W1.hat.ave = W_hat_droplast_ave(betak, Skb1, betadelta, npc)
W1.hat.ave.nopen = W_hat_droplast_ave_nopen(betak, Skb1, betadelta, npc)
W1.hat.ave.nopen = W_hat_droplast_ave_nopen(betak, Skb1, betadelta, npc)

W2.hat.nopen = W_hat_droplast_nopen(betak2, Skb2, npc)
W2.hat.ave = W_hat_droplast_ave(betak2, Skb2, betadelta2, npc)
W2.hat.ave.nopen = W_hat_droplast_ave_nopen(betak2, Skb2, betadelta2, npc)

#------------------------------------------
# CV for adaptive soft-thresholding SSL1
#------------------------------------------

W1.hat = W1.hat.nopen
W2.hat = W2.hat.nopen

betaSSL = (betadelta1 + betadelta2 
           - drop(W1.hat %*% Sk2)-drop(W2.hat %*% Sk1))/2
betaSSLk = (betak1 + betak2 
            - W1.hat %*% t(Skb2) - W2.hat %*% t(Skb1))/2

betaSSL.nopen = betaSSL
betaSSLk.nopen = betaSSLk

#------------------------------------------
# CV for adaptive soft-thresholding SSL2
#------------------------------------------

W1.hat = W1.hat.ave
W2.hat = W2.hat.ave

betaSSL = (betadelta1 + betadelta2 
           - drop(W1.hat %*% Sk2)-drop(W2.hat %*% Sk1))/2
betaSSLk = (betak1 + betak2 
            - W1.hat %*% t(Skb2) - W2.hat %*% t(Skb1))/2

betaSSL.ave = betaSSL
betaSSLk.ave = betaSSLk

#------------------------------------------
# CV for adaptive soft-thresholding SSL3
#------------------------------------------

W1.hat = W1.hat.ave.nopen
W2.hat = W2.hat.ave.nopen

betaSSL = (betadelta1 + betadelta2 
           - drop(W1.hat %*% Sk2)-drop(W2.hat %*% Sk1))/2
betaSSLk = (betak1 + betak2 
            - W1.hat %*% t(Skb2) - W2.hat %*% t(Skb1))/2

betaSSL.ave.nopen = betaSSL
betaSSLk.ave.nopen = betaSSLk