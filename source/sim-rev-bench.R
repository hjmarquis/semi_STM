.libPaths(c("~/R-3.6.1/library",.libPaths()))

library(MASS)
library(survival)
library(doParallel)
library(doRNG)

source("source/beta_delta.R")
source("source/beta_delta_perturb.R")

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
    
    h = sd(data$C)/(sum(dataS$delta1))^0.25
    KC = dnorm(as.matrix(dist(data$C/h,diag=T,upper=T)))/h
    
    betadelta1 = init.beta(dataS$delta1,data$Z,KC,
                          link = link, dlink = dlink)
    
    h = sd(data$C)/(sum(dataS$delta2))^0.25
    KC = dnorm(as.matrix(dist(data$C/h,diag=T,upper=T)))/h
    
    betadelta2 = init.beta(dataS$delta2,data$Z,KC,
                           link = link, dlink = dlink)

    
    # betak1 =  matrix(0,p, Nperturb)
    # betak2 =  matrix(0,p, Nperturb)
    # for (iperturb in 1:Nperturb) 
    # {
    #   V = rbeta(n,0.5,1.5)*4
    #   betak1[,iperturb] = 
    #     init.beta.perturb(dataS$delta1,data$Z,KC1,
    #                       V, 
    #                       init = betadelta1,
    #                       link = link, dlink = dlink)
    #   betak2[,iperturb] = 
    #     init.beta.perturb(dataS$delta2,data$Z,KC2,
    #                       V, 
    #                       init = betadelta2,
    #                       link = link, dlink = dlink)
    # }

    save(betadelta1, betadelta2, # betak1, betak2,
         file=file.name)
    
    return(NULL)
    # }
  }


stopCluster(cl)

#------------------------------------------
# One-step estimator and its permutation
#------------------------------------------

  file.name = paste("simresult/",setup,
                    version, ".rda",sep = '')   
  beta1 = beta2 =  matrix(0,p,Nrep)
  beta1s = beta2s = vector("list",Nrep)
  for(irep in 1:Nrep) 
  {
    load(paste("simresult/",version,"/",setup,
               '_rep', 
               repid[irep], ".rda",sep = ''))
    
    beta1[,irep] = betadelta1
    beta2[,irep] = betadelta2
    beta1s[[irep]] = betak1
    beta2s[[irep]] = betak2
    
  }
  save(beta1, beta2, beta1s, beta2s, 
       file = file.name)
