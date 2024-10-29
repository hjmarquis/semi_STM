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
    label.dat = sim.gen.rev(m,beta0,sigmaZ,m, misT = misT, inv.link = inv.link)$data
    tmpdat = sim.gen.rev(n,beta0,sigmaZ,m, misT = misT, inv.link = inv.link)
    data = tmpdat$data
    dataS = tmpdat$dataS
    
    #------------------------------------------
    #    beta_delta and its perturbation
    #------------------------------------------
    
    h = sd(label.dat$C[1:m])/(sum(label.dat$delta[1:m]))^0.25
    KC = dnorm(as.matrix(dist(label.dat$C[1:m]/h,diag=T,upper=T)))/h
    
    betadelta = init.beta(label.dat$delta[1:m],label.dat$Z[1:m,],KC,
                          link = link, dlink = dlink)

    betak =  matrix(0,p, Nperturb)
    for (iperturb in 1:Nperturb) 
    {
      tmpdat = sim.gen.rev(m,beta0,sigmaZ,m, misT = misT, inv.link = inv.link)$data
      betak[,iperturb] = init.beta(tmpdat$delta[1:m],tmpdat$Z[1:m,],
                                   dnorm(as.matrix(dist(tmpdat$C[1:m]/h,diag=T,upper=T)))/h,
                                   link = link, dlink = dlink)
    }
    
    betak2 =  matrix(0,p, Nperturb)
    for (iperturb in 1:Nperturb) 
    {
      tmpdat = sim.gen.rev(m,beta0,sigmaZ,m, misT = misT, inv.link = inv.link)$data
      betak2[,iperturb] = init.beta(tmpdat$delta[1:m],tmpdat$Z[1:m,],
                                   dnorm(as.matrix(dist(tmpdat$C[1:m]/h,diag=T,upper=T)))/h,
                                   link = link, dlink = dlink)
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
    
    for (iperturb in 1:Nperturb) 
    {
      tmpdat = sim.gen.rev(n,beta0,sigmaZ,m, misT = misT, inv.link = inv.link)
      data = tmpdat$data
      dataS = tmpdat$dataS
      
      betak.std = betak[,iperturb]/sqrt(sum(betak[,iperturb]^2))
      lp = drop(tmpdat$data$Z %*% betak.std)
      Skb[iperturb,1:p] = drop(Sk_sym(lp,tmpdat$data$Z,
                                      tmpdat$dataS$X1,tmpdat$dataS$delta1,
                                      tmpdat$data$C,dnorm,h1))
      Skb[iperturb,p+1:p] = drop(Sk_sym(lp,tmpdat$data$Z,
                                        tmpdat$dataS$X2,tmpdat$dataS$delta2,
                                        tmpdat$data$C,dnorm,h2))
    }
    
    Skb2 =  matrix(0, Nperturb,p*2)
    for (iperturb in 1:Nperturb) 
    {
      tmpdat = sim.gen.rev(n,beta0,sigmaZ,m, misT = misT, inv.link = inv.link)
      data = tmpdat$data
      dataS = tmpdat$dataS
      
      betak.std = betak2[,iperturb]/sqrt(sum(betak[,iperturb]^2))
      lp = drop(tmpdat$data$Z %*% betak.std)
      Skb2[iperturb,1:p] = drop(Sk_sym(lp,tmpdat$data$Z,
                                      tmpdat$dataS$X1,tmpdat$dataS$delta1,
                                      tmpdat$data$C,dnorm,h1))
      Skb2[iperturb,p+1:p] = drop(Sk_sym(lp,tmpdat$data$Z,
                                        tmpdat$dataS$X2,tmpdat$dataS$delta2,
                                        tmpdat$data$C,dnorm,h2))
    }
    
    #------------------------------------------
    #    Estimating W matrix
    #------------------------------------------
    
    W.hat.nopen = W_hat_droplast_nopen(betak, Skb, npc)
    W.hat.ave = W_hat_droplast_ave(betak, Skb, betadelta, npc)
    W.hat.ave.nopen = W_hat_droplast_ave_nopen(betak, Skb, betadelta, npc)
    
    #------------------------------------------
    # CV for adaptive soft-thresholding SSL1
    #------------------------------------------
    
    W.hat = W.hat.nopen
    
    betaSSL = betadelta - drop(W.hat %*% Sk)
    betaSSLk = betak2 - W.hat %*% t(Skb2)
    
    betaSSL.nopen = betaSSL
    betaSSLk.nopen = betaSSLk
    
    #------------------------------------------
    # CV for adaptive soft-thresholding SSL2
    #------------------------------------------

    W.hat = W.hat.ave
    
    betaSSL = betadelta - drop(W.hat %*% Sk)
    betaSSLk = betak2 - W.hat %*% t(Skb2)
    
    betaSSL.ave = betaSSL
    betaSSLk.ave = betaSSLk
    
    #------------------------------------------
    # CV for adaptive soft-thresholding SSL3
    #------------------------------------------
    
    W.hat = W.hat.ave.nopen
    
    betaSSL = betadelta - drop(W.hat %*% Sk)
    betaSSLk = betak2 - W.hat %*% t(Skb2)
    
    betaSSL.ave.nopen = betaSSL
    betaSSLk.ave.nopen = betaSSLk
    
    save(betadelta, betak, betak2, 
         W.hat.ave, betaSSL.ave,betaSSLk.ave, 
         W.hat.ave.nopen, betaSSL.ave.nopen,betaSSLk.ave.nopen, 
         W.hat.nopen, betaSSL.nopen,betaSSLk.nopen, 
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
  betaG = beta.init =  matrix(0,p,Nrep)
  betaGs = betas.init =  vector("list",Nrep)
  for(irep in 1:Nrep) 
  {
    load(paste("simresult/",version,"/",setup,
               '_rep', 
               repid[irep], ".rda",sep = ''))
    
    beta.init[,irep] = betadelta
    betas.init[[irep]] = betak2
    betaG[,irep] = betaSSL.nopen
    betaGs[[irep]] = betaSSLk.nopen
  }
  save(betaG, betaGs, beta.init,betas.init,
       file = file.name)



  file.name = paste("simresult/",setup,
                    version, "_ave.rda",sep = '')   
  betaG = beta.init =   matrix(0,p,Nrep)
  betaGs = betas.init = vector("list",Nrep)
  for(irep in 1:Nrep) 
  {
    load(paste("simresult/",version,"/",setup,
               '_rep', 
               repid[irep], ".rda",sep = ''))
    
    beta.init[,irep] = betadelta
    betas.init[[irep]] = betak
    betaG[,irep] = betaSSL.ave
    betaGs[[irep]] = betaSSLk.ave
    
  }
  save(betaG, betaGs, beta.init,betas.init,
       file = file.name)

for (isetup in 1:Nsetup)
{
  file.name = paste("simresult/",setup[isetup],
                    version, "_ave_nopen.rda",sep = '')   
  betaG = beta.init =  matrix(0,p,Nrep)
  betaGs = betas.init = vector("list",Nrep)
  for(irep in 1:Nrep) 
  {
    load(paste("simresult/",version,"/",setup[isetup],
               '_rep', 
               repid[irep], ".rda",sep = ''))
    
    beta.init[,irep] = betadelta
    betas.init[[irep]] = betak
    betaG[,irep] = betaSSL.ave.nopen
    betaGs[[irep]] = betaSSLk.ave.nopen
    
  }
  save(betaG, betaGs, beta.init,betas.init,
       file = file.name)
}