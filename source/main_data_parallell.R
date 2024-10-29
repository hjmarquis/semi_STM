# Standardize the covariates
#------------------------------------------
p = ncol(data$Z)
npc = p-1

#    beta_delta and its perturbation
#------------------------------------------

  h = sd(data$C[1:m])/(sum(data$delta[1:m]))^0.25
  KC = dnorm(as.matrix(dist(data$C[1:m]/h,diag=T,upper=T)))/h
  
  std.weight = data$IPW[1:m]/mean(data$IPW[1:m])
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
  
  betadelta = init.beta.perturb(data$delta[1:m],data$Z[1:m,],KC,
                                std.weight,
                                 init = tmpinit,
                                  maxit = 1000,
                                min.factor = 1, tol=1e-8)


#    Sk
#------------------------------------------

  # Sk = rep(0,p*2)
  beta.std = betadelta/sqrt(sum(betadelta^2))
  lp = drop(data$Z %*% beta.std)
  sdlp = sd(lp)
  h1 = sdlp/(sum(data$delta1))^0.3
  # Sk[1:p] = Sk_sym_perturb(matrix(lp),data$Z,data$X1,data$delta1,
  #                  data$C,dnorm,h1, matrix(data$IPW))
  # Sk[1:p] = Sk_sym(matrix(lp),data$Z,data$X1,data$delta1,
  #                  data$Cstar,dnorm,h1)
  h2 = sdlp/(sum(data$delta2))^0.3
  # Sk[p+1:p] = Sk_sym_perturb(matrix(lp),data$Z,data$X2,data$delta2,
  #                    data$C,dnorm,h2, matrix(data$IPW))
  # Sk[p+1:p] = Sk_sym(matrix(lp),data$Z,data$X2,data$delta2,
  #                    data$Cstar,dnorm,h2)

#    perturbation
#------------------------------------------
  
  h = sd(data$C[1:m])/(sum(data$delta[1:m]))^0.25
  KC = dnorm(as.matrix(dist(data$C[1:m]/h,diag=T,upper=T)))/h
  tmpout<- foreach(iperturb = 0:Nperturb,.packages = c("MASS","glmnet")) %dopar%
  {
    if(iperturb > 0)
    {
      V = rbeta(n,0.5,1.5)*4
   
      div = TRUE
      
      while(div)
      {
        div = FALSE
        std.weight = (data$IPW[1:m]*V[1:m])/mean(data$IPW[1:m]*V[1:m])
        tmpinit = glm(delta~log(C)+Z,data=data[1:m,],family=binomial,
                      weights = std.weight)$coef[3:(p+2)]
        fit = try(betak <-
                    init.beta.perturb(data$delta[1:m],data$Z[1:m,],KC,
                                      std.weight, 
                                      # init = betadelta,
                                      init = tmpinit,
                                      maxit = 1000,
                                      min.factor = 1, tol=1e-9))
        
        if(class(fit) == "try-error")
        {
          print(paste("Error in iteration ", iperturb))
          div = TRUE
          V[1:m] = rbeta(m,0.5,1.5)*4
        }
      } 
    }else
    {
      V = rep(1,n)
      betak = betadelta
    }
    betak.std = betak/sqrt(sum(betak^2))
    lp = data$Z %*% betak.std
    
    Skb = c(Sk_sym_perturb(lp,data$Z,data$X1,data$delta1,
                               data$Cstar,dnorm,h1,
                               matrix(V)),
            Sk_sym_perturb(lp,data$Z,data$X2,data$delta2,
                                 data$Cstar,dnorm,h2,
                                 matrix(V)))
    return(list(betak = betak, Skb=Skb))
  }
  
  Sk =  tmpout[[1]]$Skb
  
  Skb =  do.call(rbind,sapply(tmpout[-1], '[',"Skb"))
  betak =  do.call(cbind,sapply(tmpout[-1], '[',"betak"))

  rm(tmpout)
  
  #------------------------------------------
  # CV for adaptive soft-thresholding initial
  #------------------------------------------ 
  
  Nfold = 5
  Nlam = 100
  min.lam = 1e-4
  KC = dnorm(as.matrix(dist(data$C[1:m]/h,diag=T,upper=T)))/h
  std.weight = data$IPW[1:m]/mean(data$IPW[1:m])
  
  init.cv = matrix(0,p,Nfold)
  foldid = c(rep_len(1:Nfold,m)[sample(m)],
             rep_len(Nfold:1,n-m)[sample(n-m)])
  div = TRUE
  while (div) 
  {
    div = FALSE
    foldid[1:m] = rep_len(1:Nfold,m)[sample(m)]
    lp.init.cv = matrix(0,n,Nfold)
    for(ifold in 1:Nfold)
    {
      i.label = which(foldid[1:m] != ifold)
      i.train = which(foldid != ifold)
      tmpinit = glm(delta~log(C)+Z,data=data[i.label,],family=binomial,
                    weights = std.weight[i.label])$coef[3:(p+2)]
      
      fit = try(init.cv[,ifold] <- init.beta.perturb(data$delta[i.label],data$Z[i.label,],KC[i.label,i.label],
                                                     std.weight[i.label],
                                                     init = tmpinit,
                                                     maxit = 1000,
                                                     min.factor = 1, tol=1e-8))
      if(class(fit) == "try-error")
      {
        print(paste("Error in fold ", ifold))
        div = TRUE
        break
      }
      lp.init.cv[,ifold] = drop(data$Z %*% init.cv[,ifold])
    }
  }
  
  ada.factor = 1/(abs(betadelta)*apply(data$Z, 2, sd))
  lambda = c(0,exp(seq(log(min.lam),log(min(max(betadelta^2),m^(-.75)))
                       ,length.out = Nlam))[-Nlam])
  
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
                                  data$X1[i.test],data$delta1[i.test],
                                  data$Cstar[i.test]) + 
                     cv.concordance(lp.test, 
                                    data$X2[i.test],data$delta2[i.test],
                                    data$Cstar[i.test]))
    cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
  }
  
  best.ilam = which.max(apply(cv.concord,2,mean))
  lambda.best = lambda[best.ilam]
  thres.best = lambda.best * ada.factor
  betadelta.thres = sign(betadelta) * pmax(0,abs(betadelta) - thres.best)
  
  
#    Estimating W matrix
#------------------------------------------
  W.hat = W_hat_droplast_ave_nopen(betak, Skb, betadelta.thres, npc)

# CV for adaptive soft-thresholding
#------------------------------------------

  
  betaG = betadelta - drop(W.hat %*% Sk)
  ada.factor = 1/(abs(betaG)*apply(data$Z, 2, sd))
  lambda = c(0,exp(seq(log(min.lam),log(max(betaG^2)),length.out = Nlam-1)))
  
  cv.concord = matrix(0,Nfold,Nlam)
  for(ifold in 1:Nfold)
  {
    i.train = which(foldid != ifold)
    i.test = which(foldid == ifold)
    Sk.cv = rep(0,p*2)
    Sk.cv[1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                        data$X1[i.train],data$delta1[i.train],
                        data$Cstar[i.train],dnorm,h1)
    Sk.cv[p+1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                          data$X2[i.train],data$delta2[i.train],
                          data$Cstar[i.train],dnorm,h2)
    beta.cv = init.cv[,ifold] - drop (W.hat %*% Sk.cv)
    beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
    lp.test = data$Z[i.test,] %*% beta.thres.cv
    tmp.concord = (cv.concordance(lp.test, 
                                  data$X1[i.test],data$delta1[i.test],
                                  data$Cstar[i.test]) + 
                     cv.concordance(lp.test, 
                                    data$X2[i.test],data$delta2[i.test],
                                    data$Cstar[i.test]))
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

save(betaG, betaGs, beta.init,betas.init,
     betaG.thres,betaGs.thres,betaGs.thres2,
     betaG.thres.std, betaGs.thres.std,betaGs.thres2.std,
     Z.sd, 
     file = paste(out.dir,"final-",version,".rda",sep=''))

