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
Sk11 = rep(0,p*2)
Sk11[1:p] = Sk_sym(lp1[fold1],data$Z[fold1,],
                  dataS$X1[fold1],dataS$delta1[fold1],
                  data$C[fold1],dnorm,h1)
Sk11[p+1:p] = Sk_sym(lp1[fold1],data$Z[fold1,],
                    dataS$X2[fold1],dataS$delta2[fold1],
                    data$C[fold1],dnorm,h2)
Sk12 = rep(0,p*2)
Sk12[1:p] = Sk_sym(lp2[fold1],data$Z[fold1,],
                   dataS$X1[fold1],dataS$delta1[fold1],
                   data$C[fold1],dnorm,h1)
Sk12[p+1:p] = Sk_sym(lp2[fold1],data$Z[fold1,],
                     dataS$X2[fold1],dataS$delta2[fold1],
                     data$C[fold1],dnorm,h2)

Sk22 = rep(0,p*2)
Sk22[1:p] = Sk_sym(lp2[fold2],data$Z[fold2,],
                  dataS$X1[fold2],dataS$delta1[fold2],
                  data$C[fold2],dnorm,h1)
Sk22[p+1:p] = Sk_sym(lp2[fold2],data$Z[fold2,],
                    dataS$X2[fold2],dataS$delta2[fold2],
                    data$C[fold2],dnorm,h2)
Sk21 = rep(0,p*2)
Sk22[1:p] = Sk_sym(lp1[fold2],data$Z[fold2,],
                   dataS$X1[fold2],dataS$delta1[fold2],
                   data$C[fold2],dnorm,h1)
Sk22[p+1:p] = Sk_sym(lp1[fold2],data$Z[fold2,],
                     dataS$X2[fold2],dataS$delta2[fold2],
                     data$C[fold2],dnorm,h2)

#------------------------------------------
#    Perturbation 1: for Wk
#------------------------------------------

betak1 =  matrix(0,p, Nperturb)    
betak2 =  matrix(0,p, Nperturb)
Skb11 = Skb12 = matrix(0, Nperturb,p*2)
Skb21 = Skb22 = matrix(0, Nperturb,p*2)     
for (iperturb in 1:Nperturb) 
{
  V = rbeta(n,0.5,1.5)*4
  betak1[,iperturb] = init.beta.perturb(data$delta[fold1.m],data$Z[fold1.m,],KC1,
                                        V[fold1.m], 
                                        init = betadelta1,
                                        link = link, dlink = dlink)
  betak2[,iperturb] = init.beta.perturb(data$delta[fold2.m],data$Z[fold2.m,],KC2,
                                        V[fold2.m], 
                                        init = betadelta2,
                                        link = link, dlink = dlink)
  
  betak1.std = betak1[,iperturb]/sqrt(sum(betak1[,iperturb]^2))
  lp1 = (data$Z %*% betak1.std)
  betak2.std = betak2[,iperturb]/sqrt(sum(betak2[,iperturb]^2))
  lp2 = (data$Z %*% betak2.std)
  
  Skb11[iperturb,1:p] = drop(Sk_sym_perturb(lp1[fold1,,drop=F],data$Z[fold1,],
                                           dataS$X1[fold1],dataS$delta1[fold1],
                                           data$C[fold1],dnorm,h1,
                                           matrix(V[fold1])))
  Skb11[iperturb,p+1:p] = drop(Sk_sym_perturb(lp1[fold1,,drop=F],data$Z[fold1,],
                                             dataS$X2[fold1],dataS$delta2[fold1],
                                             data$C[fold1],dnorm,h2,
                                             matrix(V[fold1])))  
  Skb12[iperturb,1:p] = drop(Sk_sym_perturb(lp2[fold1,,drop=F],data$Z[fold1,],
                                            dataS$X1[fold1],dataS$delta1[fold1],
                                            data$C[fold1],dnorm,h1,
                                            matrix(V[fold1])))
  Skb12[iperturb,p+1:p] = drop(Sk_sym_perturb(lp2[fold1,,drop=F],data$Z[fold1,],
                                              dataS$X2[fold1],dataS$delta2[fold1],
                                              data$C[fold1],dnorm,h2,
                                              matrix(V[fold1])))
  
  
  Skb22[iperturb,1:p] = drop(Sk_sym_perturb(lp2[fold2,,drop=F],data$Z[fold2,],
                                           dataS$X1[fold2],dataS$delta1[fold2],
                                           data$C[fold2],dnorm,h1,
                                           matrix(V[fold2])))
  Skb22[iperturb,p+1:p] = drop(Sk_sym_perturb(lp2[fold2,,drop=F],data$Z[fold2,],
                                             dataS$X2[fold2],dataS$delta2[fold2],
                                             data$C[fold2],dnorm,h2,
                                             matrix(V[fold2])))  
  Skb21[iperturb,1:p] = drop(Sk_sym_perturb(lp1[fold2,,drop=F],data$Z[fold2,],
                                            dataS$X1[fold2],dataS$delta1[fold2],
                                            data$C[fold2],dnorm,h1,
                                            matrix(V[fold2])))
  Skb21[iperturb,p+1:p] = drop(Sk_sym_perturb(lp1[fold2,,drop=F],data$Z[fold2,],
                                              dataS$X2[fold2],dataS$delta2[fold2],
                                              data$C[fold2],dnorm,h2,
                                              matrix(V[fold2])))
}

#------------------------------------------
#    Estimating W matrix
#------------------------------------------

W1.hat.nopen = W_hat_droplast_nopen(betak1, Skb1, npc)
W1.hat.ave = W_hat_droplast_ave(betak1, Skb1, betadelta1, npc)
W1.hat.ave.nopen = W_hat_droplast_ave_nopen(betak1, Skb1, betadelta1, npc)

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