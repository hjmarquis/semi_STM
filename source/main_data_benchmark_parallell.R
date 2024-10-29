# Standardize the covariates
#------------------------------------------
p = ncol(data$Z)
npc = p-1

#    beta_delta and its perturbation
#------------------------------------------

  h1 = sd(data$C)/(sum(data$delta1))^0.25
  KC1 = dnorm(as.matrix(dist(data$C/h1,diag=T,upper=T)))/h1
  
  betadelta1 = init.beta(data$delta1,data$Z,KC1,
                         link = link, dlink = dlink)
  
  h2 = sd(data$C)/(sum(data$delta2))^0.25
  KC2 = dnorm(as.matrix(dist(data$C/h2,diag=T,upper=T)))/h2
  
  betadelta2 = init.beta(data$delta2,data$Z,KC2,
                         link = link, dlink = dlink)
  
  tmpout<- foreach(iperturb = 0:Nperturb,.packages = c("MASS","glmnet")) %dopar%
  {     
    V = rbeta(n,0.5,1.5)*4
  
    div = TRUE
  
    while(div)
    {
      div = FALSE
      fit1 = try(tmpbeta1 <-
                   init.beta.perturb(data$delta1,data$Z,KC1,
                                     V, 
                                     init = betadelta1,
                                     link = link, dlink = dlink))
      fit2 = try(tmpbeta2 <-
                   init.beta.perturb(data$delta1,data$Z,KC1,
                                     V, 
                                     init = betadelta1,
                                     link = link, dlink = dlink))
    
      if((class(fit1) == "try-error") | (class(fit2) == "try-error"))
      {
        print(paste("Error in iteration ", iperturb))
        div = TRUE
        V = rbeta(n,0.5,1.5)*4
      }
    }
    return(list(betak1 = tmpbeta1, betak2=tmpbeta2))
  }
  betak1 =  do.call(cbind,sapply(tmpout[-1], '[',"betak1"))
  betak2 =  do.call(cbind,sapply(tmpout[-1], '[',"betak2"))


save(betadelta1, betadelta2, betak1, betak2, 
     file = paste(out.dir,"final-",version,"bench.rda",sep=''))

