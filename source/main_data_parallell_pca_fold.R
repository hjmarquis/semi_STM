
#    beta_delta and its perturbation
#------------------------------------------

betadelta.fold <- foreach(ifold = 1:neg.nfold, .combine = "cbind") %dopar%
  {
    data$delta = data[[paste0("delta",ifold)]]
    data$IPW = data[[paste0("IPW",ifold)]]
    fold.pos = which(data$IPW>0)
    h = wtd.sd(data$C, data$IPW)/(sum(data$delta[fold.pos]))^0.25
    KC = dnorm(as.matrix(dist(data$C[fold.pos]/h,diag=T,upper=T)))/h
    
    tmpinit = glm(delta~log(C)+Z,data=data,family=binomial,
                  weights = IPW)$coef[3:(p+2)]
    
    if(any(is.na(tmpinit)))
    {
      data$Z = data$Z[,!is.na(tmpinit)]
      Z.sd = Z.sd[!is.na(tmpinit)]
      p = ncol(data$Z)
      npc = p-1 
      tmpinit = glm(delta~log(C)+Z,data=data,family=binomial,
                    weights = IPW)$coef[3:(p+2)]
    }
    
    start.time = Sys.time()
    betadelta = init.beta.perturb(data$delta[fold.pos],data$Z[fold.pos,],KC,
                                  data$IPW[fold.pos],
                                  init = tmpinit,
                                  maxit = 1000,
                                  min.factor = 1, tol=1e-8)
    run.time = Sys.time() - start.time
    
    return(betadelta)
  }
betadelta = apply(betadelta.fold, 1, mean)


#    Sk
#------------------------------------------

  Sk = rep(0,p*2)
  beta.std = betadelta/sqrt(sum(betadelta^2))
  lp = drop(data$Z %*% beta.std)
  sdlp = sd(lp)
  h1 = sdlp/(sum(dataS$delta1))^0.3
  Sk[1:p] = Sk_sym(matrix(lp),data$Z,dataS$X1,dataS$delta1,
                   data$C,dnorm,h1)
  h2 = sdlp/(sum(dataS$delta2))^0.3
  Sk[p+1:p] = Sk_sym(matrix(lp),data$Z,dataS$X2,dataS$delta2,
                     data$C,dnorm,h2)

#    perturbation
#------------------------------------------
  
  iperturb = 1
  tmpout<- foreach(iperturb = 1:Nperturb,.packages = c("MASS","glmnet")) %dopar%
  {
    V = rbeta(n,0.5,1.5)*4
    div = TRUE
    betak.fold = betadelta.fold
    while(div)
    {
      div = FALSE
      for(ifold in 1:neg.nfold)
      {
        data$delta = data[[paste0("delta",ifold)]]
        data$IPW = data[[paste0("IPW",ifold)]]
        fold.pos = which(data$IPW>0)
        h = wtd.sd(data$C, data$IPW)/(sum(data$delta[fold.pos]))^0.25
        KC = dnorm(as.matrix(dist(data$C[fold.pos]/h,diag=T,upper=T)))/h
        
        std.weight = (data$IPW[fold.pos]*V[fold.pos])/mean(data$IPW[fold.pos]*V[fold.pos])
        fit = try(betak <-
                    init.beta.perturb(data$delta[fold.pos],data$Z[fold.pos,],KC,
                                      std.weight, 
                                      init = betadelta.fold[,ifold],
                                      maxit = 1000,
                                      min.factor = 1, tol=1e-9))
        
        if(class(fit) == "try-error")
        {
          print(paste("Error in iteration ", iperturb))
          div = TRUE
          V = rbeta(n,0.5,1.5)*4
          break
        }
        betak.fold[,ifold] = betak
      } 
    }
  
    betak = apply(betak.fold,1,mean)
    betak.std = betak/sqrt(sum(betak^2))
    lp = data$Z %*% betak.std
    
    Skb = c(Sk_sym_perturb(lp,data$Z,dataS$X1,dataS$delta1,
                               data$C,dnorm,h1,
                               matrix(V)),
            Sk_sym_perturb(lp,data$Z,dataS$X2,dataS$delta2,
                                 data$C,dnorm,h2,
                                 matrix(V)))
    return(list(betak = betak, Skb=Skb))
  }
  
  Skb =  do.call(rbind,sapply(tmpout[-1], '[',"Skb"))
  betak =  do.call(cbind,sapply(tmpout[-1], '[',"betak"))

  rm(tmpout)
  
  
  #    Estimating W matrix
  #------------------------------------------
  W.hat = W_hat_adaPCA_ridge(betak, Skb)$W.hat
  
  #------------------------------------------
  # CV for adaptive soft-thresholding initial
  #------------------------------------------ 
  
  Nfold = 5
  Nlam = 100
  min.lam = 1e-4
  
  init.cv = matrix(0,p,Nfold)
  foldid = rep(0,n)
  pos1 = which(data$delta==1)

  div = TRUE
  while (div) 
  {
    div = FALSE
    foldid[pos1] = rep_len(1:Nfold,length(pos1))[sample(length(pos1))]
    foldid[-pos1] = rep_len(1:Nfold,n-length(pos1))[sample(n-length(pos1))]
    lp.init.cv = matrix(0,n,Nfold)
    for(ifold in 1:Nfold)
    {
      # i.label = which(foldid[1:m] != ifold)
      i.train = foldid != ifold
      
      fit = try(init.cv[,ifold] <-
                  apply( foreach(jfold = 1:neg.nfold, .combine = "cbind") %dopar%
                          {
                            data$delta = data[[paste0("delta",jfold)]]
                            data$IPW = data[[paste0("IPW",jfold)]]
                            fold.pos = which(data$IPW>0)
                            h = wtd.sd(data$C, data$IPW)/(sum(data$delta[fold.pos]))^0.25
                            data$IPW = data[[paste0("IPW",jfold)]]*i.train
                            fold.pos = which(data$IPW>0)
                            KC = dnorm(as.matrix(dist(data$C[fold.pos]/h,diag=T,upper=T)))/h
                            
                            start.time = Sys.time()
                            tmp.beta = init.beta.perturb(data$delta[fold.pos],data$Z[fold.pos,],KC,
                                                          data$IPW[fold.pos],
                                                          init = betadelta.fold[,jfold],
                                                          maxit = 1000,
                                                          min.factor = 1, tol=1e-8)
                            run.time = Sys.time() - start.time
                            
                            return(betadelta)
                          },1,mean)
                  )
      
      if(class(fit) == "try-error")
      {
        print(paste("Error in fold ", ifold))
        div = TRUE
        break
      }
      lp.init.cv[,ifold] = drop(data$Z %*% init.cv[,ifold])
    }
  }
  # 
  # ada.factor = 1/(abs(betadelta)*apply(data$Z, 2, sd))
  # lambda = c(0,exp(seq(log(min.lam),log(min(max(betadelta^2),m^(-.75)))
  #                      ,length.out = Nlam))[-Nlam])
  # 
  # #beta.cv = matrix(0,p,Nfold)
  # cv.concord = matrix(0,Nfold,Nlam)
  # for(ifold in 1:Nfold)
  # {
  #   i.train = which(foldid != ifold)
  #   i.test = which(foldid == ifold)
  #   beta.cv = init.cv[,ifold] 
  #   beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
  #   lp.test = data$Z[i.test,] %*% beta.thres.cv
  #   tmp.concord = (cv.concordance(lp.test, 
  #                                 data$X1[i.test],data$delta1[i.test],
  #                                 data$Cstar[i.test]) + 
  #                    cv.concordance(lp.test, 
  #                                   data$X2[i.test],data$delta2[i.test],
  #                                   data$Cstar[i.test]))
  #   cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
  # }
  # 
  # best.ilam = which.max(apply(cv.concord,2,mean))
  # lambda.best = lambda[best.ilam]
  # thres.best = lambda.best * ada.factor
  # betadelta.thres = sign(betadelta) * pmax(0,abs(betadelta) - thres.best)
  

# CV for adaptive soft-thresholding
#------------------------------------------

  
  betaG = betadelta - drop(W.hat %*% Sk)
  ada.factor = 1/abs(betaG)
  lambda = c(0,exp(seq(log(min.lam),log(max(betaG^2)),length.out = Nlam-1)))
  
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


# One-step estimator and its permutation
#------------------------------------------

beta.init = betadelta
betas.init = betak
betaG = (
  betadelta - drop(W.hat %*% Sk)
)
betaGs = (
  betak - W.hat %*% t(Skb)
)
thres.best = 1/abs(betaG) * lambda.best
betaG.thres = sign(betaG) * pmax(0,abs(betaG) - thres.best)
betaGs.thres = sign(betaGs) * pmax(0,abs(betaGs) - thres.best)
betaG.thres.std = (betaG.thres/
                     sqrt(sum(betaG.thres^2))*
                     sqrt(sum(betaG^2)))

betaGs.thres.std = (betaGs.thres/
                      rep(sqrt(apply(betaGs.thres^2,2,sum))/
                            sqrt(apply(betaGs^2,2,sum)),each = p))

thresGs = 1/abs(betaGs) * lambda.best
betaGs.thres2 = sign(betaGs) * pmax(0,abs(betaGs) - thresGs)
betaGs.thres2.std = (betaGs.thres2/
                       rep(sqrt(apply(betaGs.thres2^2,2,sum))/
                             sqrt(apply(betaGs^2,2,sum)),each = p))



