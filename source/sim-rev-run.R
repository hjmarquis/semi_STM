.libPaths(c("~/R-3.6.1/library",.libPaths()))

library(MASS)
library(survival)
library(doParallel)
library(doRNG)

source("source/beta_delta.R")
source("source/beta_delta_perturb.R")
source("source/Sk_sym.R")
source("source/Sk_sym_perturb.R")
source("source/W_hat_droplast_nopen.R")
source("source/W_hat_droplast_ave.R")
source("source/W_hat_droplast_ave_nopen.R")
source("source/cv_concordance.R")

np = 20
cl = makeCluster(getOption("cl.cores", np))
registerDoParallel(cl)
registerDoRNG(seed = 531)
clusterEvalQ(cl, .libPaths(c("~/R-3.6.1/library",.libPaths())))

#------------------------------------------
#    Data Generation Event Time
#------------------------------------------

nullout <- foreach(irep = 1:Nrep,.packages = c("MASS","glmnet")) %dopar%
  {
    file.name = paste("simresult/",version,"/",setup[isetup],
                      '_rep', 
                      repid[irep], ".rda",sep = '')
    if(file.exists(file.name))
      return(NULL)
    tmp = sim.gen.rev(n,beta0,sigmaZ,m, misT = misT, inv.link = inv.link)
    foldid = c(rep_len(1:Nfold,m)[sample(m)],
               rep_len(Nfold:1,n-m)[sample(n-m)])
    data = tmp$data
    dataS = tmp$dataS
    
    #------------------------------------------
    #    beta_delta and its perturbation
    #------------------------------------------
    
    h = sd(data$C[1:m])/(sum(data$delta[1:m]))^0.25
    KC = dnorm(as.matrix(dist(data$C[1:m]/h,diag=T,upper=T)))/h
    
    betadelta = init.beta(data$delta[1:m],data$Z[1:m,],KC,
                          link = link, dlink = dlink)
    
    init.cv = matrix(0,p,Nfold)
    lp.init.cv = matrix(0,n,Nfold)
    for(ifold in 1:Nfold)
    {
      i.label = which(foldid[1:m] != ifold)
      i.train = which(foldid != ifold)
      init.cv[,ifold] = init.beta(data$delta[i.label],data$Z[i.label,],KC[i.label,i.label],
                                  link = link, dlink = dlink)
      lp.init.cv[,ifold] = drop(data$Z %*% init.cv[,ifold])
    }
    
    
    
    
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
    Sk = rep(0,p*2)
    Sk[1:p] = Sk_sym(lp,data$Z,
                     dataS$X1,dataS$delta1,
                     data$C,dnorm,h1)
    Sk[p+1:p] = Sk_sym(lp,data$Z,
                       dataS$X2,dataS$delta2,
                       data$C,dnorm,h2)
    
    #------------------------------------------
    #    Sk perturbation
    #------------------------------------------
    
    
    Skb =  matrix(0, Nperturb,p*2)
    
    h1 = sdlp/(sum(dataS$delta1))^0.3
    h2 = sdlp/(sum(dataS$delta2))^0.3
    
    betak =  matrix(0,p, Nperturb)
    for (iperturb in 1:Nperturb) 
    {
      V = rbeta(n,0.5,1.5)*4
      betak[,iperturb] = 
        init.beta.perturb(data$delta[1:m],data$Z[1:m,],KC,
                          V[1:m], 
                          init = betadelta,
                          link = link, dlink = dlink)
      betak.std = betak[,iperturb]/sqrt(sum(betak[,iperturb]^2))
      lp = data$Z %*% betak.std
      Skb[iperturb,1:p] = c(Sk_sym_perturb(lp,data$Z,
                                           dataS$X1,dataS$delta1,
                                           data$C,dnorm,h1,
                                           matrix(V)))
      Skb[iperturb,p+1:p] = c(Sk_sym_perturb(lp,data$Z,
                                             dataS$X2,dataS$delta2,
                                             data$C,dnorm,h2,
                                             matrix(V)))
    }
    
    #------------------------------------------
    # CV for adaptive soft-thresholding initial
    #------------------------------------------ 
    ada.factor = 1/(abs(betadelta)*apply(data$Z, 2, sd))
    lambda = c(0,exp(seq(log(min.lam),log(max(betadelta^2)),length.out = Nlam))[-Nlam])
    
    #beta.cv = matrix(0,p,Nfold)
    cv.concord = matrix(0,Nfold,Nlam)
    for(ifold in 1:Nfold)
    {
      i.train = which(foldid != ifold)
      i.test = which(foldid == ifold)
      beta.cv = init.cv[,ifold] 
      beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
      lp.test = data$Z[i.test,] %*% beta.thres.cv
      tmp.concord = (cv.concordance(lp.test, 
                                    dataS$X1[i.test],dataS$delta1[i.test],
                                    data$C[i.test]) + 
                       cv.concordance(lp.test, 
                                      dataS$X2[i.test],dataS$delta2[i.test],
                                      data$C[i.test]))
      cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
    }
    
    best.ilam = which.max(apply(cv.concord,2,mean))
    lambda.best = lambda[best.ilam]
    thres.best = lambda.best * ada.factor
    betadelta.thres = sign(betadelta) * pmax(0,abs(betadelta) - thres.best)
    
    #------------------------------------------
    #    Estimating W matrix
    #------------------------------------------
    
    W.hat.nopen = W_hat_droplast_nopen(betak, Skb, npc)
    W.hat.ave = W_hat_droplast_ave(betak, Skb, betadelta.thres, npc)
    W.hat.ave.nopen = W_hat_droplast_ave_nopen(betak, Skb, betadelta.thres, npc)
    
    #------------------------------------------
    # CV for adaptive soft-thresholding SSL1
    #------------------------------------------
    
    W.hat = W.hat.nopen
    
    betaSSL = betadelta - drop(W.hat %*% Sk)
    ada.factor = 1/(abs(betaSSL)*apply(data$Z, 2, sd))
    lambda = c(0,exp(seq(log(min.lam),log(max(betaSSL^2)),length.out = Nlam))[-Nlam])
    
    #beta.cv = matrix(0,p,Nfold)
    cv.concord = matrix(0,Nfold,Nlam)
    for(ifold in 1:Nfold)
    {
      i.train = which(foldid != ifold)
      i.test = which(foldid == ifold)
      Sk.cv = rep(0,p*2)
      Sk.cv[1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                          dataS$X1[i.train],dataS$delta1[i.train],
                          data$C[i.train],dnorm,h1)
      Sk.cv[p+1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                            dataS$X2[i.train],dataS$delta2[i.train],
                            data$C[i.train],dnorm,h2)
      beta.cv = init.cv[,ifold] - drop (W.hat %*% Sk.cv)
      beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
      lp.test = data$Z[i.test,] %*% beta.thres.cv
      tmp.concord = (cv.concordance(lp.test, 
                                    dataS$X1[i.test],dataS$delta1[i.test],
                                    data$C[i.test]) + 
                       cv.concordance(lp.test, 
                                      dataS$X2[i.test],dataS$delta2[i.test],
                                      data$C[i.test]))
      cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
    }
    
    best.ilam = which.max(apply(cv.concord,2,mean))
    lambda.best = lambda[best.ilam]
    thres.best = lambda.best * ada.factor
    
    betaSSLk = betak - W.hat %*% t(Skb)
    
    betaSSL.nopen = betaSSL
    betaSSLk.nopen = betaSSLk
    thres.best.nopen = thres.best
    
    #------------------------------------------
    # CV for adaptive soft-thresholding SSL2
    #------------------------------------------

    W.hat = W.hat.ave
    
    betaSSL = betadelta - drop(W.hat %*% Sk)
    ada.factor = 1/(abs(betaSSL)*apply(data$Z, 2, sd))
    lambda = c(0,exp(seq(log(min.lam),log(max(betaSSL^2)),length.out = Nlam))[-Nlam])
    
    #beta.cv = matrix(0,p,Nfold)
    cv.concord = matrix(0,Nfold,Nlam)
    for(ifold in 1:Nfold)
    {
      i.train = which(foldid != ifold)
      i.test = which(foldid == ifold)
      Sk.cv = rep(0,p*2)
      Sk.cv[1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                          dataS$X1[i.train],dataS$delta1[i.train],
                          data$C[i.train],dnorm,h1)
      Sk.cv[p+1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                            dataS$X2[i.train],dataS$delta2[i.train],
                            data$C[i.train],dnorm,h2)
      beta.cv = init.cv[,ifold] - drop (W.hat %*% Sk.cv)
      beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
      lp.test = data$Z[i.test,] %*% beta.thres.cv
      tmp.concord = (cv.concordance(lp.test, 
                                    dataS$X1[i.test],dataS$delta1[i.test],
                                    data$C[i.test]) + 
                       cv.concordance(lp.test, 
                                      dataS$X2[i.test],dataS$delta2[i.test],
                                      data$C[i.test]))
      cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
    }
    
    best.ilam = which.max(apply(cv.concord,2,mean))
    lambda.best = lambda[best.ilam]
    thres.best = lambda.best * ada.factor
    
    betaSSLk = betak - W.hat %*% t(Skb)
    
    
    betaSSL.ave = betaSSL
    betaSSLk.ave = betaSSLk
    thres.best.ave = thres.best
    
    #------------------------------------------
    # CV for adaptive soft-thresholding SSL3
    #------------------------------------------
    
    W.hat = W.hat.ave.nopen
    
    betaSSL = betadelta - drop(W.hat %*% Sk)
    ada.factor = 1/(abs(betaSSL)*apply(data$Z, 2, sd))
    lambda = c(0,exp(seq(log(min.lam),log(max(betaSSL^2)),length.out = Nlam))[-Nlam])
    
    #beta.cv = matrix(0,p,Nfold)
    cv.concord = matrix(0,Nfold,Nlam)
    for(ifold in 1:Nfold)
    {
      i.train = which(foldid != ifold)
      i.test = which(foldid == ifold)
      Sk.cv = rep(0,p*2)
      Sk.cv[1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                          dataS$X1[i.train],dataS$delta1[i.train],
                          data$C[i.train],dnorm,h1)
      Sk.cv[p+1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                            dataS$X2[i.train],dataS$delta2[i.train],
                            data$C[i.train],dnorm,h2)
      beta.cv = init.cv[,ifold] - drop (W.hat %*% Sk.cv)
      beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
      lp.test = data$Z[i.test,] %*% beta.thres.cv
      tmp.concord = (cv.concordance(lp.test, 
                                    dataS$X1[i.test],dataS$delta1[i.test],
                                    data$C[i.test]) + 
                       cv.concordance(lp.test, 
                                      dataS$X2[i.test],dataS$delta2[i.test],
                                      data$C[i.test]))
      cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
    }
    
    best.ilam = which.max(apply(cv.concord,2,mean))
    lambda.best = lambda[best.ilam]
    thres.best = lambda.best * ada.factor
    
    betaSSLk = betak - W.hat %*% t(Skb)
    
    
    betaSSL.ave.nopen = betaSSL
    betaSSLk.ave.nopen = betaSSLk
    thres.best.ave.nopen = thres.best
    
    save(betadelta, betak, betadelta.thres, 
         W.hat.ave, betaSSL.ave,betaSSLk.ave,thres.best.ave, 
         W.hat.ave.nopen, betaSSL.ave.nopen,betaSSLk.ave.nopen,thres.best.ave.nopen, 
         W.hat.nopen, betaSSL.nopen,betaSSLk.nopen,thres.best.nopen, 
         file=file.name)
    
    return(NULL)
    # }
  }


stopCluster(cl)

#------------------------------------------
# One-step estimator and its permutation
#------------------------------------------

  file.name = paste("simresult/",setup,
                    version, "_nopen.rda",sep = '')   
  betaG = beta.init = betaG.thres = betaG.thres.std =  matrix(0,p,Nrep)
  betaGs = betas.init = betaGs.thres = betaGs.thres.std = vector("list",Nrep)
  for(irep in 1:Nrep) 
  {
    load(paste("simresult/",version,"/",setup,
               '_rep', 
               repid[irep], ".rda",sep = ''))
    
    beta.init[,irep] = betadelta
    betas.init[[irep]] = betak
    betaG[,irep] = betaSSL.nopen
    betaGs[[irep]] = betaSSLk.nopen
    betaG.thres[,irep] = sign(betaG[,irep]) * pmax(0,abs(betaG[,irep]) - thres.best.nopen)
    betaGs.thres[[irep]] = sign(betaGs[[irep]]) * pmax(0,abs(betaGs[[irep]]) - thres.best.nopen)
    betaG.thres.std[,irep] = (betaG.thres[,irep]/
                                sqrt(sum(betaG.thres[,irep]^2))*
                                sqrt(sum(betaG[,irep]^2)))
    
    betaGs.thres.std[[irep]] = (betaGs.thres[[irep]]/
                                  rep(sqrt(apply(betaGs.thres[[irep]]^2,2,sum))/
                                        sqrt(apply(betaGs[[irep]]^2,2,sum)),each = p))
    
  }
  save(betaG, betaGs, beta.init,betas.init,
       betaG.thres,betaGs.thres,
       betaG.thres.std, betaGs.thres.std,
       file = file.name)



  file.name = paste("simresult/",setup,
                    version, "_ave.rda",sep = '')   
  betaG = beta.init = betaG.thres = betaG.thres.std =  matrix(0,p,Nrep)
  betaGs = betas.init = betaGs.thres = betaGs.thres.std = vector("list",Nrep)
  for(irep in 1:Nrep) 
  {
    load(paste("simresult/",version,"/",setup,
               '_rep', 
               repid[irep], ".rda",sep = ''))
    
    beta.init[,irep] = betadelta
    betas.init[[irep]] = betak
    betaG[,irep] = betaSSL.ave
    betaGs[[irep]] = betaSSLk.ave
    betaG.thres[,irep] = sign(betaG[,irep]) * pmax(0,abs(betaG[,irep]) - thres.best.ave)
    betaGs.thres[[irep]] = sign(betaGs[[irep]]) * pmax(0,abs(betaGs[[irep]]) - thres.best.ave)
    betaG.thres.std[,irep] = (betaG.thres[,irep]/
                                sqrt(sum(betaG.thres[,irep]^2))*
                                sqrt(sum(betaG[,irep]^2)))
    
    betaGs.thres.std[[irep]] = (betaGs.thres[[irep]]/
                                  rep(sqrt(apply(betaGs.thres[[irep]]^2,2,sum))/
                                        sqrt(apply(betaGs[[irep]]^2,2,sum)),each = p))
    
  }
  save(betaG, betaGs, beta.init,betas.init,
       betaG.thres,betaGs.thres,
       betaG.thres.std, betaGs.thres.std,
       file = file.name)




for (isetup in 1:Nsetup)
{
  file.name = paste("simresult/",setup[isetup],
                    version, "_ave_nopen.rda",sep = '')   
  betaG = beta.init = betaG.thres = betaG.thres.std =  matrix(0,p,Nrep)
  betaGs = betas.init = betaGs.thres = betaGs.thres.std = vector("list",Nrep)
  for(irep in 1:Nrep) 
  {
    load(paste("simresult/",version,"/",setup[isetup],
               '_rep', 
               repid[irep], ".rda",sep = ''))
    
    beta.init[,irep] = betadelta
    betas.init[[irep]] = betak
    betaG[,irep] = betaSSL.ave.nopen
    betaGs[[irep]] = betaSSLk.ave.nopen
    betaG.thres[,irep] = sign(betaG[,irep]) * pmax(0,abs(betaG[,irep]) - thres.best.ave.nopen)
    betaGs.thres[[irep]] = sign(betaGs[[irep]]) * pmax(0,abs(betaGs[[irep]]) - thres.best.ave.nopen)
    betaG.thres.std[,irep] = (betaG.thres[,irep]/
                                sqrt(sum(betaG.thres[,irep]^2))*
                                sqrt(sum(betaG[,irep]^2)))
    
    betaGs.thres.std[[irep]] = (betaGs.thres[[irep]]/
                                  rep(sqrt(apply(betaGs.thres[[irep]]^2,2,sum))/
                                        sqrt(apply(betaGs[[irep]]^2,2,sum)),each = p))
    
  }
  save(betaG, betaGs, beta.init,betas.init,
       betaG.thres,betaGs.thres,
       betaG.thres.std, betaGs.thres.std,
       file = file.name)
}