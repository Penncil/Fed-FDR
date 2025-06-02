library(glmnet)
library(mvtnorm)
library(ncvreg)


expit<-function(x){
  1/(1+exp(-x))
}

####data generate function######
data.generate<-function(cov.type=1,
                        betashift=TRUE,
                        K=15,
                        n=200,
                        p=500,
                        Gamma=-1,
                        Corr=rep(0.3,K),
                        supp=50,
                        signal=4){
  site <- rep(1:K, each = n)
  #theta <- c(rnorm(supp,0,signal*sqrt(supp*log(p)/n)), rep(0,p-supp-1))
  si<-sample(c(-1,1), supp, replace = TRUE)
  #theta <- c(rep(signal*sqrt(log(p)/n),supp)*si, rep(0,p-supp-1))
  #theta <- c(rep(signal*sqrt(supp*log(p)/n),supp)*sample(c(-1,1), supp, replace = TRUE), rep(0,p-supp-1))
  X<-Treat<-NULL
  if(cov.type==1){
    sigmax<-diag(rep(1,p-1))
    X<-rmvnorm(n*K, rep(0,p-1), sigmax)
  }else if(cov.type==2){
    for (i in 1:K) {
      sigmax<-toeplitz(Corr[i]^(0:(p-2)))
      X1<-rmvnorm(n, rep(0,p-1), sigmax)
      X  <-rbind(X, X1)
      if(betashift){
      theta <- c(runif(supp,0.5,4)*rep(signal*sqrt(log(p)/n),supp)*si, rep(0,p-supp-1))
      }else{
        theta <- c(rep(signal*sqrt(log(p)/n),supp)*si, rep(0,p-supp-1))  
      }
      meanTreat <- expit(Gamma[i]+X1%*%theta)
      Treat<-c(Treat, rbinom(n, 1, meanTreat))
    }
  }else if(cov.type==3){
    for (i in 1:K) {
      sigmax  <- matrix(corr, p-1, p-1) + diag(1-corr, p-1)
      X  <-rbind(X, rmvnorm(n, rep(0,p-1), sigmax))
    }
  }
  
  #meanTreat <- expit(Gamma+X%*%theta)
  #Treat <- rbinom(n*K, 1, meanTreat)
  #Y <- X%*%theta+rnorm(n*K)
  
  theta.true<-c(gamma,theta)
  
  return(list(#Yall=Y,
              Treatall=Treat, 
              Xall=X, site=site,
              theta.true=theta.true))    
}



#######refined debiased Lasso####
REF_DS_inf <- function(x, y, family, lasso_est,delta) {
  nn <- length(y)
  X <- cbind(rep(1, nrow(x)), x)
  p<-ncol(X)
  if(p != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "binomial") {
    mu <- as.vector(exp(X%*%lasso_est)/(1+exp(X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%((mu*(1-mu))*X)/nn
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
  } else {
    stop("Input family is not supported.")
  }
  
  theta_inv <-tryCatch({
     solve(neg_ddloglik_glmnet)
  },error=function(e){ tryCatch({
    solve(neg_ddloglik_glmnet+delta*diag(rep(1,p)))
  },error=function(e){
    cat("delta is too small \n")
  })
   }
  )
  b_hat_inv <- as.vector(lasso_est - theta_inv%*%neg_dloglik_glmnet)
  se_inv <- sqrt(diag(theta_inv))/sqrt(nn)
  pval_inv <- 2*pnorm(abs(b_hat_inv/se_inv), lower.tail=F)
  
  return(list(est=b_hat_inv, se=se_inv, pvalue=pval_inv, theta=theta_inv))
}

#### node wised lasso
ORIG_DS_inf <- function(x, y, family, lasso_est, nfold=5, n_lambda=100, lambda_ratio=0.005) {
  nn <- length(y)
  pp <- ncol(x)
  X <- cbind(rep(1, nrow(x)), x)
  if(ncol(X) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "binomial") {
    mu <- as.vector(1/(1+exp(-X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu*(1-mu))%*%X/nn
    C_glmnet <- sqrt(diag(mu*(1-mu))/nn)%*%X
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
    C_glmnet <- sqrt(diag(mu)/nn)%*%X
  } else {
    stop("Input family is not supported.")
  }
  
  theta_glmnet <- diag(pp+1)
  tau_glmnet <- rep(NA, pp+1)
  for(j in 1:(pp+1)) { # for: nodewise lasso
    current_x <- sqrt(nn)*C_glmnet[,-j]
    current_y <- sqrt(nn)*as.vector(C_glmnet[,j])
    lam_max <- max(abs(t(current_x)%*%current_y)/nn)
    lam_min <- lam_max*lambda_ratio
    lam_seq <- exp(seq(from=log(lam_max), to=log(lam_min), length.out=n_lambda))
    gamma_j_glmnet <- cv.glmnet(x=current_x, y=current_y,
                                family="gaussian", alpha=1, standardize=F, intercept=F,
                                nfolds=nfold, lambda=lam_seq)
    gamma_j_glmnet <- as.vector(glmnet(x=sqrt(n)*C_glmnet[,-j], y=sqrt(n)*as.vector(C_glmnet[,j]),
                                       family="gaussian", alpha=1, standardize=F, intercept=F,
                                       lambda=gamma_j_glmnet$lambda.min)$beta)
    theta_glmnet[j,-j] <- (-1)*t(gamma_j_glmnet)
    tau_glmnet[j] <- as.numeric(neg_ddloglik_glmnet[j,j]-neg_ddloglik_glmnet[j,-j]%*%gamma_j_glmnet)
  } # end for: nodewise lasso
  theta_glmnet <- diag(1/tau_glmnet)%*%theta_glmnet 
  
  b_hat_nw <- as.vector(lasso_est - theta_glmnet%*%neg_dloglik_glmnet)
  se_nw <- sqrt(diag(theta_glmnet%*%neg_ddloglik_glmnet%*%t(theta_glmnet)))/sqrt(nn)
  pval_nw <- 2*pnorm(abs(b_hat_nw/se_nw), lower.tail=F)
  
  return(list(est=b_hat_nw, se=se_nw, pvalue=pval_nw, theta=theta_glmnet))
}


####

ORIG_DS_inf_simulation<- function(x, y, family, lasso_est) {
  nn <- length(y)
  pp <- ncol(x)
  X <- cbind(rep(1, nrow(x)), x)
  if(ncol(X) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "binomial") {
    mu <- as.vector(1/(1+exp(-X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%(mu*(1-mu)*X)/nn
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
    C_glmnet <- sqrt(diag(mu)/nn)%*%X
  } else {
    stop("Input family is not supported.")
  }
  
 
  tau_glmnet <- diag(neg_ddloglik_glmnet)
  
  theta_glmnet <- 1/tau_glmnet
  
  b_hat_nw <- as.vector(lasso_est - theta_glmnet*neg_dloglik_glmnet)
  se_nw <- theta_glmnet/nn
  pval_nw <- 2*pnorm(abs(b_hat_nw/sqrt(se_nw)), lower.tail=F)
  #se_nw1<-neg_ddloglik_glmnet%*%(neg_ddloglik_glmnet)
  
  return(list(est=b_hat_nw, se=se_nw, pvalue=pval_nw))
}

####


FDR_control<-function(beta1, beta2, q){
  M<- sign(beta1*beta2)*(abs(beta1)+abs(beta2)) 
  tau.seq<-abs(M[which(M!=0)])
  if(!any(M<0)){
    tau.hat<-0
  }else{ 
    tau.seq<-sort(tau.seq)
    for (tau in c(0,tau.seq)) {
      fdp<- sum(I(M<=(-tau)))/max(sum(I(M>=tau)),1)
      if(fdp<=q){
        tau.hat<-tau
        break
      } 
      tau.hat<-tau
    }
  }
  p<-length(M)
  index<-which(M>tau.hat)
  S<-rep(0,p)
  S[index]<-1
  return(S)
}


MFDR<-function(Theta,q){
  S<-NULL
  p<-dim(Theta)[2]
  K<-dim(Theta)[1]
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      beta1<-Theta[i,]
      beta2<-Theta[j,]
      S<-rbind(S,FDR_control(beta1, beta2,q))
    }
  }
  
  a<-rowSums(S)
  tau.seq<-colSums(S / (a*I(a>=1)+1*I(a<1)))/(K*(K-1)/2)
  index<-order(tau.seq)
  I.seq<-tau.seq[index]
  index1<-which(cumsum(I.seq)>q)
  supp<-index[index1]
  b.FDR<-rep(0,p)
  b.FDR[supp]=1
  
  return(list(FDR=b.FDR))
}





BH<-function(pvalue,q){
p<-length(pvalue)
plevel<-(1:p)*q/p

index<-order(pvalue)
p.seq<-pvalue[index]   

i=1
while (i<p) {
  if(p.seq[i]>plevel[i]){
    index1<-i-1
    break
  }
  i<-i+1
} 
supp<-index[1:index1]
b.FDR<-rep(0,p)
b.FDR[supp]=1

return(list(FDR=b.FDR)) 
}




Tpr<-function(thetahat,theta){
  tp<-sum(thetahat*theta)
  tn<-sum(thetahat*(1-theta))  
  fdp<-tn/max(sum(thetahat),1)
  power<-tp/sum(theta)
  
  return(c(fdp,power))
}


local.glm<-function(Treatall,
                    Xall,
                    site,
                    delta=1e-5){
  
  siteId<-unique(site)
  K <- length(siteId)
  nsite <- table(site)
  d<-dim(Xall)[2]
  
  Theta_local<-matrix(0, nrow = K, ncol = d)
  Lasso.es<-se<-matrix(0, nrow = K, ncol = d)
  #b1<-numeric(d+1)
  #varest<-matrix(0,d+1,d+1)
  for(i in 1:K){
    index<-which(site==siteId[i])
    Xlocal <- Xall[index,]
    Treatlocal <- Treatall[index]
    #Ylocal<-Yall[index]
    
    fit0 <- cv.glmnet(Xlocal, Treatlocal, nfolds = 5,family = "binomial")
    #lambda<-(fit0$lambda.min+fit0$lambda.1se)/2
    lambda<-fit0$lambda.min
    fit1 <- glmnet(Xlocal, Treatlocal, family = "binomial",lambda = lambda)
    Theta <- as.vector(coef(fit1))
    ###meta
    
    fit1 <- ORIG_DS_inf_simulation(x=Xlocal,y=Treatlocal, family="binomial", 
                                   Theta)
    Lasso.es[i,]<-fit1$est[-1]
    se[i,]<-1/fit1$se[-1]
    
  }
  # varest<-solve(varest)
  #estmeta<-varest%*%b1
  
  varest<-1/colSums(se)
  estmeta<-colSums(Lasso.es*se)*varest
  pvalue<-2*pnorm(abs(estmeta/sqrt(varest)),lower.tail = FALSE)
  
  return(list(pvalue=pvalue,
              delasso=Lasso.es
              #Beta=Beta_local#, Tau=Tau
  ))
}

local.glm.scad<-function(Treatall,
                         Xall,
                         site,
                         lambda_scale=1,
                         delta=1e-5){
  
  siteId<-unique(site)
  K <- length(siteId)
  nsite <- table(site)
  d<-dim(Xall)[2]
  
  Theta_local<-matrix(0, nrow = K, ncol = d)
  Lasso.es<-se<-support_all<-matrix(0, nrow = K, ncol = d)
  #b1<-numeric(d+1)
  #varest<-matrix(0,d+1,d+1)
  for(i in 1:K){
    index<-which(site==siteId[i])
    Xlocal <- Xall[index,]
    Treatlocal <- Treatall[index]
    #ridge fit
    
    fit0 <- cv.ncvreg(Xlocal, Treatlocal, nfolds = 5, family = "binomial", penalty="SCAD")
    #lambda<-(fit0$lambda.min+fit0$lambda.1se)/2
    lambda<-lambda_scale*fit0$lambda.min
    fit1 <- ncvreg(Xlocal, Treatlocal, family = "binomial", lambda = lambda, penalty="SCAD")
    Theta <- as.vector(coef(fit1))
    support_all[i,]<-I(Theta[-1]!=0)
    ###meta
  }
  # varest<-solve(varest)
  #estmeta<-varest%*%b1
  
  
  for(i in 1:K){
    index<-which(site==siteId[i])
    Xlocal <- Xall[index,]
    Treatlocal <- Treatall[index]
    
    if(K>2){
      support<-which(colSums(support_all[-i,])!=0)
    }else {
      support<-which(support_all[-i,]!=0)
    }
    X<-Xlocal[,support]
    fit0 <- cv.glmnet(X, Treatlocal,nfolds = 5, family = "binomial")
    lambda<-fit0$lambda.min
    fit0 <- glmnet(X, Treatlocal, family = "binomial",lambda = lambda)
    Theta <- as.vector(coef(fit0))
    fit0<-REF_DS_inf(X, Treatlocal, family="binomial", Theta, delta)
    Theta_local[i,support]<-as.vector(fit0$est/sqrt(fit0$se))[-1]
  }
  
  
  return(list(Theta=Theta_local
              #Beta=Beta_local#, Tau=Tau
  ))
}




local.glm.lasso<-function(Treatall,
                          Xall,
                          site,
                          lambda_scale=1,
                          weight=TRUE,
                          delta=1e-5){

  siteId<-unique(site)
  K <- length(siteId)
  nsite <- table(site)
  d<-dim(Xall)[2]

  Theta_local<-matrix(0, nrow = K, ncol = d)
  Lasso.es<-se<-support_all<-matrix(0, nrow = K, ncol = d)
  #b1<-numeric(d+1)
  #varest<-matrix(0,d+1,d+1)
  for(i in 1:K){
    index<-which(site==siteId[i])
    Xlocal <- Xall[index,]
    Treatlocal <- Treatall[index]
   if(weight==TRUE){
    ridge_fit <- glmnet(Xlocal ,  Treatlocal, family = "binomial", alpha = 0)
    beta_ridge <- coef(ridge_fit, s = 0.01)[-1]
    weights_ridge <- 1 / abs(beta_ridge)
    weights_ridge[is.infinite(weights_ridge)] <- max(weights_ridge[is.finite(weights_ridge)])

    fit0 <- cv.glmnet(Xlocal, Treatlocal, nfolds = 5, family = "binomial", penalty.factor = weights_ridge)
    #lambda<-lambda_scale*(fit0$lambda.min+fit0$lambda.1se)/2
    lambda<-lambda_scale*fit0$lambda.min
    #lambda<-lambda_scale*fit0$lambda.1se
    fit1 <- glmnet(Xlocal, Treatlocal, family = "binomial", lambda = lambda, penalty.factor = weights_ridge)
   }else{
    fit0 <- cv.glmnet(Xlocal, Treatlocal, nfolds = 5, family = "binomial")
    #lambda<-lambda_scale*(fit0$lambda.min+fit0$lambda.1se)/2
    lambda<-lambda_scale*fit0$lambda.min
    #lambda<-lambda_scale*fit0$lambda.1se
    fit1 <- glmnet(Xlocal, Treatlocal, family = "binomial", lambda = lambda)
   }
    Theta <- as.vector(coef(fit1))
    support_all[i,]<-I(Theta[-1]!=0)
  }
  #varest<-solve(varest)
  #estmeta<-varest%*%b1


  for(i in 1:K){
    index<-which(site==siteId[i])
    Xlocal <- Xall[index,]
    Treatlocal <- Treatall[index]
    if(K>2){
    support<-which(colSums(support_all[-i,])!=0)
    }else {
      support<-which(support_all[-i,]!=0)
    }
    X<-Xlocal[,support]
    fit0 <- cv.glmnet(X, Treatlocal,nfolds = 5, family = "binomial")
    lambda<-fit0$lambda.min
    fit0 <- glmnet(X, Treatlocal, family = "binomial",lambda = lambda)
    Theta <- as.vector(coef(fit0))
    fit0<-REF_DS_inf(X, Treatlocal, family="binomial", Theta, delta)
    Theta_local[i,support]<-as.vector(fit0$est/sqrt(fit0$se))[-1]
  }


  return(list(Theta=Theta_local
              #Beta=Beta_local#, Tau=Tau
  ))
}




#####BHq method

BHq<-function(Treatall,
              Xall){
    fit0 <- cv.glmnet(Xall, Treatall, nfolds = 5,family = "binomial")
    lambda<-fit0$lambda.min
    fit0 <- glmnet(Xall, Treatall, family = "binomial",lambda = lambda)
    lasso_est<- as.vector(coef(fit0))
    
    fit0 <- ORIG_DS_inf_simulation(x=Xall,y=Treatall, family="binomial",lasso_est)
    
    pvalue<-fit0$pvalue[-1]
  return(list(pvalue=pvalue
  ))
}











